"""
track1_evolutionary.py
Track 1: Evolutionary-Guided Variant Design for Maize D1

ENHANCED with Contrastive Learning:
- Positive examples: Heat-adapted species (learn FROM)
- Negative examples: Cold-adapted species (avoid traits)
- Enrichment scoring: High in positive, low in negative
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple
from pathlib import Path
from collections import Counter
import json
import sys

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

from analysis_engine import (
    REFERENCE_SEQUENCES,
    get_sequence,
    get_temperature,
    ConservationAnalysis,
    MultipleSequenceAlignment
)


class EvolutionaryGuidedDesigner:
    """
    Track 1: Design variants based on evolutionary analysis
    
    ENHANCED Strategy:
    1. Categorize sequences as positive/negative/target
    2. Find residues ENRICHED in heat-adapted (positive)
    3. Find residues DEPLETED in cold-adapted (negative)
    4. Calculate enrichment score: (% in positive) - (% in negative)
    5. Generate maize variants with high-enrichment mutations
    """
    
    def __init__(self, sequences: Dict = None, temperatures: Dict = None):
        """
        Args:
            sequences: Dict mapping species name -> sequence data (uses REFERENCE_SEQUENCES if None)
            temperatures: Dict mapping species name -> temperature (legacy support)
        """
        self.sequences = sequences or REFERENCE_SEQUENCES
        
        # Extract just sequences for compatibility
        self.sequence_strings = {}
        for species, data in self.sequences.items():
            if isinstance(data, dict):
                self.sequence_strings[species] = data['sequence']
            else:
                self.sequence_strings[species] = data
        
        # Get maize sequence
        maize_data = self.sequences['maize']
        self.maize_sequence = maize_data['sequence'] if isinstance(maize_data, dict) else maize_data
        
        # Protected regions (DO NOT MUTATE)
        self.protected_regions = [
            (130, 150, "QB binding region"),
            (160, 170, "TyrZ catalytic region"),
            (210, 230, "QB site proper"),
            (332, 344, "Mn4CaO5 cluster ligands"),
        ]
        
        # Categorize sequences
        self.positive_species = self._get_species_by_label('positive')
        self.negative_species = self._get_species_by_label('negative')
        self.target_species = self._get_species_by_label('target')
        
        # Results storage
        self.alignment = None
        self.conservation_scores = None
        self.candidate_positions = None
    
    def _get_species_by_label(self, label: str) -> List[str]:
        """Get list of species with a specific label"""
        species_list = []
        for species, data in self.sequences.items():
            if isinstance(data, dict) and data.get('label') == label:
                species_list.append(species)
        return species_list
    
    def run_complete_analysis(self, 
                            min_enrichment: float = 0.5,
                            max_conservation: float = 0.7,
                            n_variants: int = 20,
                            use_contrastive: bool = True) -> Dict:
        """
        Run complete Track 1 analysis pipeline
        
        Args:
            min_enrichment: Minimum enrichment score (positive - negative prevalence)
            max_conservation: Maximum conservation score (lower = more flexible)
            n_variants: Number of variants to generate
            use_contrastive: Use contrastive analysis (True) or simple correlation (False)
        
        Returns:
            Dict with analysis results and generated variants
        """
        
        print("\n" + "="*70)
        print("🧬 TRACK 1: EVOLUTIONARY-GUIDED VARIANT DESIGN")
        print("   ENHANCED with Contrastive Learning")
        print("="*70)
        print(f"Target: Zea mays D1 protein")
        print(f"Reference species: {len(self.sequences)}")
        
        # Show categorization
        print(f"\n📊 Sequence Classification:")
        print(f"   🎯 Target (template):        {len(self.target_species)} - {', '.join(self.target_species)}")
        print(f"   ✅ Positive (heat-adapted):  {len(self.positive_species)} - {', '.join(self.positive_species)}")
        print(f"   ❌ Negative (cold-adapted):  {len(self.negative_species)} - {', '.join(self.negative_species)}")
        print("="*70 + "\n")
        
        # Step 1: Multiple sequence alignment
        print("📊 Step 1/4: Multiple sequence alignment...")
        self.alignment = MultipleSequenceAlignment.align_multiple(self.sequence_strings)
        print(f"   ✅ Aligned {len(self.alignment)} sequences")
        
        # Step 2: Conservation analysis
        print("\n📊 Step 2/4: Conservation analysis...")
        self.conservation_scores = ConservationAnalysis.conservation_score(self.alignment)
        n_flexible = np.sum(self.conservation_scores < max_conservation)
        print(f"   ✅ Found {n_flexible} flexible positions (conservation < {max_conservation})")
        
        # Step 3: Identify candidate positions (CONTRASTIVE)
        print("\n📊 Step 3/4: Contrastive analysis (heat vs. cold adaptation)...")
        
        if use_contrastive and self.positive_species and self.negative_species:
            self.candidate_positions = self._identify_candidate_positions_contrastive(
                min_enrichment=min_enrichment,
                max_conservation=max_conservation
            )
        else:
            print("   ⚠️  Insufficient labeled data for contrastive analysis")
            print("   Falling back to simple comparison...")
            self.candidate_positions = self._identify_candidate_positions_simple(
                max_conservation=max_conservation
            )
        
        print(f"   ✅ Identified {len(self.candidate_positions)} candidate positions")
        
        if self.candidate_positions:
            print(f"\n   🏆 Top 5 candidates:")
            for i, pos in enumerate(self.candidate_positions[:5], 1):
                score_label = 'enrichment_score' if 'enrichment_score' in pos else 'score'
                print(f"      {i}. Position {pos['position']}: "
                      f"{pos['from_target']}→{pos['to_positive']} "
                      f"(Score={pos[score_label]:.2f})")
        
        # Step 4: Generate variants
        print(f"\n📊 Step 4/4: Generating {n_variants} maize variants...")
        variants = self._generate_variants(n_variants)
        print(f"   ✅ Generated {len(variants)} variants")
        
        # Compile results
        results = {
            'track': 'evolutionary_guided',
            'method': 'contrastive_learning' if use_contrastive else 'simple_comparison',
            'classification': {
                'target': self.target_species,
                'positive': self.positive_species,
                'negative': self.negative_species,
            },
            'analysis': {
                'n_sequences': len(self.sequences),
                'n_flexible_positions': int(n_flexible),
                'n_candidate_positions': len(self.candidate_positions),
            },
            'candidate_positions': self.candidate_positions,
            'variants': variants,
            'parameters': {
                'min_enrichment': min_enrichment,
                'max_conservation': max_conservation,
                'use_contrastive': use_contrastive,
            }
        }
        
        return results
    
    def _is_protected(self, position: int) -> Tuple[bool, str]:
        """Check if position is in a protected region"""
        for start, end, region_name in self.protected_regions:
            if start <= position <= end:
                return True, region_name
        return False, ""
    
    def _identify_candidate_positions_contrastive(self, 
                                                  min_enrichment: float,
                                                  max_conservation: float) -> List[Dict]:
        """
        Identify positions using contrastive analysis
        
        Enrichment Score = (Prevalence in POSITIVE) - (Prevalence in NEGATIVE)
        
        High score = mutation is enriched in heat-adapted, depleted in cold-adapted
        """
        
        candidates = []
        target_seq = self.maize_sequence
        
        for pos in range(len(target_seq)):
            # Skip protected regions
            is_protected, region_name = self._is_protected(pos + 1)
            if is_protected:
                continue
            
            # Check conservation
            if pos < len(self.conservation_scores):
                conservation = self.conservation_scores[pos]
                if conservation > max_conservation:
                    continue
            else:
                continue
            
            target_aa = target_seq[pos]
            
            # Get residues in positive examples (heat-adapted)
            positive_residues = []
            for sp in self.positive_species:
                seq = self.sequence_strings[sp]
                if pos < len(seq) and seq[pos] != '-':
                    positive_residues.append(seq[pos])
            
            # Get residues in negative examples (cold-adapted)
            negative_residues = []
            for sp in self.negative_species:
                seq = self.sequence_strings[sp]
                if pos < len(seq) and seq[pos] != '-':
                    negative_residues.append(seq[pos])
            
            if not positive_residues:
                continue
            
            # Find consensus in positive examples
            positive_counts = Counter(positive_residues)
            positive_aa, pos_count = positive_counts.most_common(1)[0]
            
            # Only suggest if different from target
            if positive_aa == target_aa:
                continue
            
            # Calculate prevalence scores
            positive_prevalence = pos_count / len(positive_residues)
            negative_prevalence = negative_residues.count(positive_aa) / max(len(negative_residues), 1) if negative_residues else 0
            
            # Enrichment = high in positive, low in negative
            enrichment_score = positive_prevalence - negative_prevalence
            
            # Only keep if enriched in positives
            if enrichment_score < min_enrichment:
                continue
            
            # Bonus if negatives kept the target AA (validates we should change it)
            negative_has_target = negative_residues.count(target_aa) / max(len(negative_residues), 1) if negative_residues else 0
            if negative_has_target > 0.5:
                enrichment_score += 0.2
            
            candidates.append({
                'position': pos + 1,
                'from_target': target_aa,
                'to_positive': positive_aa,
                'enrichment_score': enrichment_score,
                'positive_prevalence': positive_prevalence,
                'negative_prevalence': negative_prevalence,
                'negative_has_target': negative_has_target,
                'conservation': float(conservation),
                'mutation': f"{target_aa}{pos+1}{positive_aa}",
                'found_in_positive': [sp for sp in self.positive_species if self.sequence_strings[sp][pos] == positive_aa],
                'found_in_negative': [sp for sp in self.negative_species if self.sequence_strings[sp][pos] == positive_aa] if negative_residues else [],
                'rationale': self._generate_contrastive_rationale(
                    pos+1, target_aa, positive_aa,
                    positive_prevalence, negative_prevalence,
                    [sp for sp in self.positive_species if self.sequence_strings[sp][pos] == positive_aa],
                    [sp for sp in self.negative_species if self.sequence_strings[sp][pos] == positive_aa] if negative_residues else []
                )
            })
        
        # Sort by enrichment score
        candidates.sort(key=lambda x: x['enrichment_score'], reverse=True)
        
        return candidates
    
    def _identify_candidate_positions_simple(self, max_conservation: float) -> List[Dict]:
        """
        Simple comparison (fallback if insufficient labeled data)
        """
        candidates = []
        target_seq = self.maize_sequence
        
        # Just compare to all non-target sequences
        other_species = [sp for sp in self.sequences.keys() if sp not in self.target_species]
        
        for pos in range(len(target_seq)):
            is_protected, _ = self._is_protected(pos + 1)
            if is_protected:
                continue
            
            if pos < len(self.conservation_scores):
                conservation = self.conservation_scores[pos]
                if conservation > max_conservation:
                    continue
            else:
                continue
            
            target_aa = target_seq[pos]
            
            # Find most common non-target residue
            other_residues = []
            for sp in other_species:
                seq = self.sequence_strings[sp]
                if pos < len(seq) and seq[pos] != '-':
                    other_residues.append(seq[pos])
            
            if not other_residues:
                continue
            
            counts = Counter(other_residues)
            consensus_aa, count = counts.most_common(1)[0]
            
            if consensus_aa != target_aa:
                prevalence = count / len(other_residues)
                
                candidates.append({
                    'position': pos + 1,
                    'from_target': target_aa,
                    'to_positive': consensus_aa,
                    'score': prevalence,
                    'conservation': float(conservation),
                    'mutation': f"{target_aa}{pos+1}{consensus_aa}",
                    'rationale': f"Found in {count}/{len(other_residues)} other species"
                })
        
        candidates.sort(key=lambda x: x['score'], reverse=True)
        return candidates
    
    def _generate_contrastive_rationale(self, position: int, from_aa: str, to_aa: str,
                                       pos_prev: float, neg_prev: float,
                                       pos_species: List, neg_species: List) -> str:
        """Generate rationale for contrastive mutation"""
        
        parts = []
        
        parts.append(f"Position {position} flexible, {from_aa}→{to_aa}")
        
        if pos_species:
            parts.append(f"found in {pos_prev*100:.0f}% of heat-adapted ({', '.join(pos_species)})")
        
        if neg_species:
            parts.append(f"but only {neg_prev*100:.0f}% of cold-adapted ({', '.join(neg_species)})")
        else:
            parts.append(f"absent in cold-adapted species")
        
        enrichment = pos_prev - neg_prev
        parts.append(f"enrichment={enrichment:.2f}")
        
        return "; ".join(parts)
    
    def _generate_variants(self, n_variants: int) -> List[Dict]:
        """
        Generate maize variants by applying candidate mutations
        """
        variants = []
        
        # Generate single-mutation variants
        for i, candidate in enumerate(self.candidate_positions[:n_variants]):
            variant_seq = list(self.maize_sequence)
            pos_idx = candidate['position'] - 1
            variant_seq[pos_idx] = candidate['to_positive']
            variant_seq = ''.join(variant_seq)
            
            identity = sum(1 for a, b in zip(self.maize_sequence, variant_seq) 
                          if a == b) / len(self.maize_sequence)
            
            score_key = 'enrichment_score' if 'enrichment_score' in candidate else 'score'
            
            variants.append({
                'variant_id': f"TRACK1_EVOL_{i+1:03d}",
                'sequence': variant_seq,
                'mutations': [candidate['mutation']],
                'n_mutations': 1,
                'identity_to_wt': identity,
                'design_rationale': candidate['rationale'],
                'evidence_score': candidate[score_key],
                'source_track': 'evolutionary_guided'
            })
        
        # Generate double-mutation variants
        if len(self.candidate_positions) >= 2:
            double_count = 0
            for i in range(min(5, len(self.candidate_positions)-1)):
                for j in range(i+1, min(i+4, len(self.candidate_positions))):
                    if double_count >= 5:
                        break
                    
                    cand1 = self.candidate_positions[i]
                    cand2 = self.candidate_positions[j]
                    
                    # Ensure mutations are >15 residues apart
                    if abs(cand1['position'] - cand2['position']) < 15:
                        continue
                    
                    variant_seq = list(self.maize_sequence)
                    variant_seq[cand1['position']-1] = cand1['to_positive']
                    variant_seq[cand2['position']-1] = cand2['to_positive']
                    variant_seq = ''.join(variant_seq)
                    
                    identity = sum(1 for a, b in zip(self.maize_sequence, variant_seq) 
                                  if a == b) / len(self.maize_sequence)
                    
                    score_key = 'enrichment_score' if 'enrichment_score' in cand1 else 'score'
                    avg_score = (cand1[score_key] + cand2[score_key]) / 2
                    
                    variants.append({
                        'variant_id': f"TRACK1_EVOL_DBL_{double_count+1:02d}",
                        'sequence': variant_seq,
                        'mutations': [cand1['mutation'], cand2['mutation']],
                        'n_mutations': 2,
                        'identity_to_wt': identity,
                        'design_rationale': f"Combined: {cand1['mutation']} + {cand2['mutation']}",
                        'evidence_score': avg_score,
                        'source_track': 'evolutionary_guided'
                    })
                    double_count += 1
        
        return variants
    
    def export_results(self, output_dir: Path, results: Dict):
        """Export Track 1 results to files"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Export variants as FASTA
        fasta_file = output_dir / "track1_evolutionary_variants.fasta"
        with open(fasta_file, 'w') as f:
            for variant in results['variants']:
                mutations_str = "+".join(variant['mutations'])
                f.write(f">{variant['variant_id']} | {mutations_str} | "
                       f"Score={variant['evidence_score']:.3f}\n")
                seq = variant['sequence']
                for i in range(0, len(seq), 60):
                    f.write(f"{seq[i:i+60]}\n")
        
        print(f"\n💾 FASTA saved: {fasta_file}")
        
        # Export candidate positions as CSV
        if results['candidate_positions']:
            csv_file = output_dir / "track1_candidate_positions.csv"
            df = pd.DataFrame(results['candidate_positions'])
            df.to_csv(csv_file, index=False)
            print(f"💾 Candidates saved: {csv_file}")
        
        # Export classification report
        report_file = output_dir / "track1_contrastive_report.md"
        with open(report_file, 'w') as f:
            f.write("# Track 1: Evolutionary-Guided Design - Contrastive Analysis Report\n\n")
            f.write("## Sequence Classification\n\n")
            f.write(f"- **Target (template)**: {', '.join(results['classification']['target'])}\n")
            f.write(f"- **Positive (heat-adapted)**: {', '.join(results['classification']['positive'])}\n")
            f.write(f"- **Negative (cold-adapted)**: {', '.join(results['classification']['negative'])}\n\n")
            
            f.write("## Top Mutations (Enriched in Heat-Adapted)\n\n")
            f.write("| Rank | Mutation | Enrichment | In Heat | In Cold | Rationale |\n")
            f.write("|------|----------|------------|---------|---------|----------|\n")
            
            for i, cand in enumerate(results['candidate_positions'][:15], 1):
                score_key = 'enrichment_score' if 'enrichment_score' in cand else 'score'
                pos_prev = cand.get('positive_prevalence', 0)
                neg_prev = cand.get('negative_prevalence', 0)
                
                f.write(f"| {i} | {cand['mutation']} | {cand[score_key]:.2f} | "
                       f"{pos_prev*100:.0f}% | {neg_prev*100:.0f}% | "
                       f"{cand['rationale'][:60]}... |\n")
        
        print(f"💾 Report saved: {report_file}")
        
        # Export full analysis as JSON
        json_file = output_dir / "track1_analysis_results.json"
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"💾 Analysis saved: {json_file}")


# Command-line interface
if __name__ == "__main__":
    from datetime import datetime
    
    print("="*70)
    print("🌽 MAIZE D1 TRACK 1: EVOLUTIONARY-GUIDED DESIGN (ENHANCED)")
    print("="*70)
    
    # Initialize with reference sequences
    designer = EvolutionaryGuidedDesigner()
    
    # Run complete analysis with contrastive learning
    results = designer.run_complete_analysis(
        min_enrichment=0.5,       # Minimum enrichment score
        max_conservation=0.7,     # Maximum conservation
        n_variants=20,
        use_contrastive=True      # Use contrastive analysis
    )
    
    # Export results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(f"../data/results/track1_{timestamp}")
    designer.export_results(output_dir, results)
    
    print("\n" + "="*70)
    print("✅ TRACK 1 COMPLETE!")
    print("="*70)
    print(f"📊 Generated {len(results['variants'])} variants")
    print(f"📁 Results saved to: {output_dir}")
    print("="*70)