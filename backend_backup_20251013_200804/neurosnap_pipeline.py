"""
neurosnap_pipeline.py
Streamlined D1 protein engineering pipeline using NeuroSnap's full suite
"""

import asyncio
import os
from typing import List, Dict, Optional
from pathlib import Path
from datetime import datetime
import json

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False
    print("⚠️  BioPython not installed. Install with: pip install biopython")

from neurosnap_client import NeurosnapClient


class StreamlinedD1Pipeline:
    """
    Simplified D1 protein engineering pipeline
    
    4 Stages (vs 7 in old pipeline):
    1. Structure Prediction (AlphaFold2)
    2. Variant Generation (RFdiffusion OR comparative genomics)
    3. Thermostability Scoring (TemStaPro)
    4. Mutation Validation (ThermoMPNN)
    """
    
    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.getenv('NEUROSNAP_API_KEY')
        if not self.api_key:
            raise ValueError("NEUROSNAP_API_KEY not found in environment")
        
        self.client = None
        self.results_dir = Path("results_v2_neurosnap")
        self.results_dir.mkdir(exist_ok=True)
    
    async def run_complete_pipeline(self,
                                   wild_type_sequence: str,
                                   target_temperature: int = 50,
                                   num_variants: int = 50,
                                   use_ai_design: bool = False) -> Dict:
        """
        Run the complete streamlined pipeline
        
        Args:
            wild_type_sequence: Starting D1 sequence
            target_temperature: Target operating temperature (°C)
            num_variants: Number of variants to generate
            use_ai_design: Use RFdiffusion (True) or comparative genomics (False)
        
        Returns:
            Complete results dictionary with top variants
        """
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        print("\n" + "="*70)
        print("🧬 STREAMLINED D1 PROTEIN ENGINEERING PIPELINE")
        print("="*70)
        print(f"📅 Timestamp: {timestamp}")
        print(f"🌡️  Target Temperature: {target_temperature}°C")
        print(f"🔢 Variants to Generate: {num_variants}")
        print(f"🤖 AI Design: {'RFdiffusion' if use_ai_design else 'Comparative Genomics'}")
        print("="*70 + "\n")
        
        async with NeurosnapClient(self.api_key) as client:
            self.client = client
            
            # STAGE 1: Predict Structure
            print("📊 STAGE 1/4: Structure Prediction (AlphaFold2)")
            print("-" * 70)
            print("⏳ This may take 5-15 minutes...")
            structure_result = await self._stage1_predict_structure(wild_type_sequence)
            
            if structure_result.get('error'):
                print(f"⚠️  Structure prediction failed, continuing with sequence-only methods")
                structure_pdb = None
            else:
                structure_pdb = structure_result.get('pdb_content')
                confidence = structure_result.get('confidence', 0)
                print(f"✅ Structure predicted with confidence: {confidence:.1f}/100")
            
            print()
            
            # STAGE 2: Generate Variants
            print("🧬 STAGE 2/4: Variant Generation")
            print("-" * 70)
            
            if use_ai_design and structure_pdb:
                variants = await self._stage2_generate_variants_ai(
                    structure_pdb, 
                    num_variants
                )
            else:
                variants = self._stage2_generate_variants_simple(
                    wild_type_sequence,
                    num_variants
                )
            
            print(f"✅ Generated {len(variants)} variants")
            print()
            
            # STAGE 3: Score Thermostability
            print("🌡️  STAGE 3/4: Thermostability Scoring (TemStaPro)")
            print("-" * 70)
            scored_variants = await self._stage3_score_thermostability(
                variants,
                target_temperature
            )
            print(f"✅ Scored {len(scored_variants)} variants")
            print()
            
            # STAGE 4: Validate Mutations
            print("✅ STAGE 4/4: Mutation Validation (ThermoMPNN)")
            print("-" * 70)
            
            if structure_pdb:
                final_variants = await self._stage4_validate_mutations(
                    wild_type_sequence,
                    scored_variants[:10],  # Only top 10 for validation
                    structure_pdb
                )
            else:
                print("⚠️  Skipping ThermoMPNN (no structure available)")
                final_variants = scored_variants
            
            print(f"✅ Validated {len(final_variants)} final variants")
            print()
            
            # Save Results
            results = self._compile_results(
                wild_type_sequence,
                final_variants,
                target_temperature,
                timestamp
            )
            
            self._save_results(results, timestamp)
            
            return results
    
    async def _stage1_predict_structure(self, sequence: str) -> Dict:
        """Stage 1: Predict 3D structure using AlphaFold2"""
        return await self.client.predict_structure_alphafold(
            sequence,
            job_name="D1_wild_type_structure"
        )
    
    async def _stage2_generate_variants_ai(self,
                                          template_pdb: str,
                                          num_designs: int) -> List[Dict]:
        """Stage 2A: Generate variants using RFdiffusion"""
        
        design_regions = [
            (50, 70),
            (120, 140),
            (200, 220),
            (280, 300)
        ]
        
        print(f"🎯 Targeting regions: {design_regions}")
        
        results = await self.client.generate_variants_rfdiffusion(
            template_pdb,
            design_regions,
            num_designs,
            job_name="D1_rfdiffusion_designs"
        )
        
        variants = []
        for i, result in enumerate(results):
            variants.append({
                'variant_id': f"RFD_{i+1:03d}",
                'sequence': result['sequence'],
                'design_score': result.get('confidence', 0),
                'method': 'rfdiffusion'
            })
        
        return variants
    
    def _stage2_generate_variants_simple(self,
                                        wild_type_sequence: str,
                                        num_variants: int) -> List[Dict]:
        """Stage 2B: Simple mutation generation (fallback)"""
        print("📚 Using simple mutation strategy...")
        
        import random
        
        # Common thermostabilizing mutations
        mutations = {
            'G': ['A', 'V'],  # Glycine to more rigid
            'S': ['T', 'A'],  # Serine substitutions
            'A': ['V', 'I'],  # Alanine to larger hydrophobic
            'N': ['D', 'Q'],  # Asparagine alternatives
            'Q': ['E', 'L'],  # Glutamine alternatives
        }
        
        variants = []
        used_mutations = set()
        
        for i in range(min(num_variants, 100)):
            # Pick a random position
            pos = random.randint(10, len(wild_type_sequence) - 10)
            from_aa = wild_type_sequence[pos]
            
            if from_aa in mutations:
                to_aa = random.choice(mutations[from_aa])
                mutation_key = f"{from_aa}{pos}{to_aa}"
                
                if mutation_key not in used_mutations:
                    used_mutations.add(mutation_key)
                    
                    # Apply mutation
                    variant_seq = list(wild_type_sequence)
                    variant_seq[pos] = to_aa
                    variant_seq = ''.join(variant_seq)
                    
                    variants.append({
                        'variant_id': f"SIMPLE_{i+1:03d}",
                        'sequence': variant_seq,
                        'mutation': mutation_key,
                        'method': 'simple_rules'
                    })
        
        return variants
    
    async def _stage3_score_thermostability(self,
                                           variants: List[Dict],
                                           target_temp: int) -> List[Dict]:
        """Stage 3: Score all variants with TemStaPro"""
        
        sequences = [v['sequence'] for v in variants]
        
        print(f"⏳ Submitting {len(sequences)} sequences to TemStaPro...")
        
        thermo_results = await self.client.batch_predict_temstapro(
            sequences,
            temperature=target_temp
        )
        
        scored_variants = []
        for variant, thermo in zip(variants, thermo_results):
            prediction = thermo.get('prediction', {})
            variant['is_stable'] = prediction.get('is_stable', False)
            variant['thermo_confidence'] = prediction.get('confidence', 0)
            variant['predicted_tm'] = prediction.get('predicted_tm', 0)
            variant['thermo_score'] = prediction.get('confidence', 0) if prediction.get('is_stable') else 0
            scored_variants.append(variant)
        
        scored_variants.sort(key=lambda x: x['thermo_score'], reverse=True)
        
        print("\n🏆 Top 5 Variants by Thermostability:")
        for i, v in enumerate(scored_variants[:5], 1):
            print(f"   {i}. {v['variant_id']}: Tm = {v['predicted_tm']:.1f}°C "
                  f"(confidence: {v['thermo_confidence']:.2f})")
        
        return scored_variants
    
    async def _stage4_validate_mutations(self,
                                        wild_type: str,
                                        variants: List[Dict],
                                        structure_pdb: str) -> List[Dict]:
        """Stage 4: Validate mutations with ThermoMPNN"""
        
        print(f"⏳ Validating top {len(variants)} variants with ThermoMPNN...")
        
        validated = []
        for variant in variants:
            mutations = self._identify_mutations(wild_type, variant['sequence'])
            
            if not mutations:
                continue
            
            mutation_effects = await self.client.predict_mutation_effects_thermompnn(
                structure_pdb,
                mutations,
                job_name=f"validate_{variant['variant_id']}"
            )
            
            if not mutation_effects or 'error' in mutation_effects:
                continue
            
            total_ddg = sum(m['predicted_ddg'] for m in mutation_effects)
            avg_confidence = sum(m['confidence'] for m in mutation_effects) / len(mutation_effects)
            
            variant['mutations_validated'] = mutation_effects
            variant['total_ddg'] = total_ddg
            variant['mutation_confidence'] = avg_confidence
            variant['is_stabilizing'] = total_ddg < 0
            
            variant['final_score'] = (
                variant['thermo_score'] * 0.5 +
                (1 - abs(total_ddg) / 10) * 0.3 +
                avg_confidence * 0.2
            )
            
            validated.append(variant)
        
        validated.sort(key=lambda x: x['final_score'], reverse=True)
        
        print("\n🏆 Top 5 Validated Variants:")
        for i, v in enumerate(validated[:5], 1):
            stabilizing = "✅" if v['is_stabilizing'] else "⚠️"
            print(f"   {i}. {v['variant_id']}: {stabilizing} ΔΔG = {v['total_ddg']:.2f} kcal/mol, "
                  f"Final Score = {v['final_score']:.3f}")
        
        return validated
    
    def _identify_mutations(self, wt_seq: str, variant_seq: str) -> List[Dict]:
        """Identify all mutations between wild-type and variant"""
        mutations = []
        for i, (wt_aa, var_aa) in enumerate(zip(wt_seq, variant_seq)):
            if wt_aa != var_aa:
                mutations.append({
                    'position': i + 1,
                    'from': wt_aa,
                    'to': var_aa
                })
        return mutations
    
    def _compile_results(self,
                        wild_type: str,
                        variants: List[Dict],
                        target_temp: int,
                        timestamp: str) -> Dict:
        """Compile final results"""
        return {
            'pipeline_version': 'v2_neurosnap_streamlined',
            'timestamp': timestamp,
            'wild_type_sequence': wild_type,
            'target_temperature': target_temp,
            'total_variants_generated': len(variants),
            'top_variants': variants[:10],
            'summary': {
                'best_variant': variants[0] if variants else None,
                'avg_tm_improvement': sum(v.get('predicted_tm', 0) - target_temp 
                                         for v in variants) / len(variants) if variants else 0,
                'num_stabilizing': sum(1 for v in variants if v.get('is_stabilizing', False))
            }
        }
    
    def _save_results(self, results: Dict, timestamp: str):
        """Save results to files"""
        
        json_file = self.results_dir / f"d1_neurosnap_{timestamp}_results.json"
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"\n💾 Results saved to: {json_file}")
        
        if HAS_BIOPYTHON:
            fasta_file = self.results_dir / f"d1_neurosnap_{timestamp}_synthesis.fasta"
            self._export_fasta(results['top_variants'], fasta_file)
            print(f"💾 FASTA saved to: {fasta_file}")
    
    def _export_fasta(self, variants: List[Dict], output_file: Path):
        """Export variants to FASTA format"""
        if not HAS_BIOPYTHON:
            return
        
        records = []
        for v in variants:
            header = (f"{v['variant_id']} | "
                     f"Tm={v.get('predicted_tm', 0):.1f}C | "
                     f"Score={v.get('final_score', v.get('thermo_score', 0)):.3f}")
            
            record = SeqRecord(
                Seq(v['sequence']),
                id=v['variant_id'],
                description=header
            )
            records.append(record)
        
        SeqIO.write(records, output_file, "fasta")


async def main():
    """Main execution function"""
    
    # Example D1 sequence (replace with your actual sequence)
    D1_WILD_TYPE = "MTAILERRESESLWGRRCNWITSTENLYI"  # Shortened for example
    
    pipeline = StreamlinedD1Pipeline()
    
    results = await pipeline.run_complete_pipeline(
        wild_type_sequence=D1_WILD_TYPE,
        target_temperature=50,
        num_variants=20,
        use_ai_design=False
    )
    
    print("\n" + "="*70)
    print("🎉 PIPELINE COMPLETED!")
    print("="*70)


if __name__ == "__main__":
    asyncio.run(main())