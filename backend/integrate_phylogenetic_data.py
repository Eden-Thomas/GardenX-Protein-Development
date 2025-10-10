# backend/integrate_phylogenetic_data.py
"""
Integrates phylogenetic analysis results into the master pipeline
Generates variant recommendations and prepares for synthesis
"""

import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Tuple
import sys

# Add backend to path
sys.path.append(str(Path(__file__).parent))

try:
    from analysis_engine import REFERENCE_SEQUENCES, calculate_conservation_score
except ImportError:
    print("⚠️  Could not import from analysis_engine - will use standalone mode")
    REFERENCE_SEQUENCES = {}


class PhylogeneticVariantGenerator:
    """
    Generates corn D1 variants based on phylogenetic analysis
    """
    
    def __init__(self, mutations_file: str = "data/phylogenetic_mutations.csv"):
        """Initialize with phylogenetic mutations data"""
        
        self.mutations_file = mutations_file
        self.mutations_df = None
        self.corn_reference = None
        
    def load_data(self):
        """Load phylogenetic mutations and reference sequence"""
        
        print("📂 Loading phylogenetic analysis results...")
        
        # Load mutations
        if Path(self.mutations_file).exists():
            self.mutations_df = pd.read_csv(self.mutations_file)
            print(f"✓ Loaded {len(self.mutations_df)} mutations")
        else:
            print(f"❌ Mutations file not found: {self.mutations_file}")
            print("   Run phylogenetic_analysis.py first!")
            return False
        
        # Load corn reference sequence
        if REFERENCE_SEQUENCES:
            # Try to find corn sequence
            for name, seq in REFERENCE_SEQUENCES.items():
                if 'corn' in name.lower() or 'zea' in name.lower() or 'mays' in name.lower():
                    self.corn_reference = seq
                    print(f"✓ Found corn reference: {name}")
                    break
        
        if not self.corn_reference:
            print("⚠️  No corn reference found - using first sequence as reference")
            if REFERENCE_SEQUENCES:
                self.corn_reference = list(REFERENCE_SEQUENCES.values())[0]
        
        return True
    
    def rank_mutations(self) -> pd.DataFrame:
        """
        Rank mutations by priority for engineering
        
        Scoring factors:
        1. Phylogenetic level (closer = safer)
        2. Conservation score
        3. Frequency across species
        4. Known stabilizing effect
        """
        
        print("\n🎯 Ranking mutations by engineering priority...")
        
        if self.mutations_df is None:
            return None
        
        # Add scoring columns
        ranked = self.mutations_df.copy()
        
        # Score 1: Phylogenetic level (closer = higher score)
        # Level 1 (corn variants) = 100 points
        # Level 2 (close C4) = 80 points
        # Level 3 = 60 points, etc.
        level_scores = {1: 100, 2: 80, 3: 60, 4: 40, 5: 20, 6: 10}
        ranked['level_score'] = ranked['level'].map(level_scores)
        
        # Score 2: Conservation (already in data)
        ranked['conservation_score'] = ranked.get('conservation_score', 0.5) * 50
        
        # Score 3: Frequency (count how many species have this mutation)
        mutation_counts = ranked.groupby(['position', 'mutant_aa']).size()
        ranked['frequency_score'] = ranked.apply(
            lambda row: mutation_counts.get((row['position'], row['mutant_aa']), 1) * 10,
            axis=1
        )
        
        # Total score
        ranked['total_score'] = (
            ranked['level_score'] + 
            ranked['conservation_score'] + 
            ranked['frequency_score']
        )
        
        # Sort by total score
        ranked = ranked.sort_values('total_score', ascending=False)
        
        return ranked
    
    def generate_variant_combinations(self, 
                                     top_n: int = 10,
                                     max_mutations_per_variant: int = 5) -> List[Dict]:
        """
        Generate variant combinations for testing
        
        Strategy:
        - Variant 1: Single best mutation (safest)
        - Variant 2-5: Top 2-5 mutations (incremental)
        - Variant 6-10: Combinations from different levels
        """
        
        print(f"\n🧬 Generating variant combinations...")
        
        ranked = self.rank_mutations()
        if ranked is None or ranked.empty:
            print("❌ No mutations to generate variants from")
            return []
        
        variants = []
        
        # Top mutations by level
        level1_muts = ranked[ranked['level'] == 1].head(5)
        level2_muts = ranked[ranked['level'] == 2].head(5)
        
        # Variant 1: Single best mutation (Level 1)
        if not level1_muts.empty:
            best = level1_muts.iloc[0]
            variants.append({
                "variant_id": "V1_Single_L1",
                "description": "Single best mutation from corn variants",
                "mutations": [
                    f"{best['reference_aa']}{best['position']}{best['mutant_aa']}"
                ],
                "expected_improvement": "Minimal risk, natural corn variation",
                "priority": 1,
                "source_level": 1
            })
        
        # Variant 2: Top 2 mutations (Level 1)
        if len(level1_muts) >= 2:
            top2 = level1_muts.head(2)
            variants.append({
                "variant_id": "V2_Double_L1",
                "description": "Top 2 mutations from corn variants",
                "mutations": [
                    f"{row['reference_aa']}{row['position']}{row['mutant_aa']}"
                    for _, row in top2.iterrows()
                ],
                "expected_improvement": "Low risk, natural combinations",
                "priority": 2,
                "source_level": 1
            })
        
        # Variant 3: Top 3 mutations (Level 1)
        if len(level1_muts) >= 3:
            top3 = level1_muts.head(3)
            variants.append({
                "variant_id": "V3_Triple_L1",
                "description": "Top 3 mutations from corn variants",
                "mutations": [
                    f"{row['reference_aa']}{row['position']}{row['mutant_aa']}"
                    for _, row in top3.iterrows()
                ],
                "expected_improvement": "Moderate improvement, validated in nature",
                "priority": 3,
                "source_level": 1
            })
        
        # Variant 4: Best from Level 2 (Sorghum, etc)
        if not level2_muts.empty:
            best_l2 = level2_muts.iloc[0]
            variants.append({
                "variant_id": "V4_Single_L2",
                "description": "Best mutation from close C4 relatives",
                "mutations": [
                    f"{best_l2['reference_aa']}{best_l2['position']}{best_l2['mutant_aa']}"
                ],
                "expected_improvement": "Proven in heat-adapted crops",
                "priority": 4,
                "source_level": 2
            })
        
        # Variant 5: Combo Level 1 + Level 2
        if not level1_muts.empty and not level2_muts.empty:
            best_l1 = level1_muts.iloc[0]
            best_l2 = level2_muts.iloc[0]
            variants.append({
                "variant_id": "V5_Combo_L1_L2",
                "description": "Best from corn + best from C4 relatives",
                "mutations": [
                    f"{best_l1['reference_aa']}{best_l1['position']}{best_l1['mutant_aa']}",
                    f"{best_l2['reference_aa']}{best_l2['position']}{best_l2['mutant_aa']}"
                ],
                "expected_improvement": "Synergistic effect expected",
                "priority": 5,
                "source_level": "1+2"
            })
        
        # Variant 6-10: Top ranked regardless of level
        top_10_overall = ranked.head(10)
        for i, group_size in enumerate([2, 3, 4, 5, 5], start=6):
            if len(variants) >= 10:
                break
            
            muts = top_10_overall.head(group_size)
            variants.append({
                "variant_id": f"V{i}_Top{group_size}",
                "description": f"Top {group_size} mutations overall",
                "mutations": [
                    f"{row['reference_aa']}{row['position']}{row['mutant_aa']}"
                    for _, row in muts.iterrows()
                ],
                "expected_improvement": f"Highest scoring {group_size}-mutation set",
                "priority": i,
                "source_level": "mixed"
            })
        
        print(f"✅ Generated {len(variants)} variant combinations")
        return variants
    
    def apply_mutations_to_sequence(self, 
                                    mutations: List[str], 
                                    reference_seq: str) -> str:
        """
        Apply mutations to reference sequence
        
        Args:
            mutations: List like ["T45A", "S127P"]
            reference_seq: Original sequence
        
        Returns:
            Modified sequence
        """
        
        seq = list(reference_seq)
        
        for mutation in mutations:
            # Parse mutation like "T45A"
            ref_aa = mutation[0]
            position = int(mutation[1:-1])
            new_aa = mutation[-1]
            
            # Apply mutation (position is 1-indexed)
            if position <= len(seq):
                if seq[position-1] == ref_aa:
                    seq[position-1] = new_aa
                else:
                    print(f"⚠️  Mismatch at position {position}: expected {ref_aa}, found {seq[position-1]}")
        
        return ''.join(seq)
    
    def generate_synthesis_files(self, variants: List[Dict], output_dir: str = "results_phylogenetic"):
        """
        Generate FASTA files ready for gene synthesis
        """
        
        print(f"\n📄 Generating synthesis-ready files in {output_dir}/...")
        
        Path(output_dir).mkdir(exist_ok=True)
        
        # Protein sequences
        protein_fasta = Path(output_dir) / "corn_d1_variants_protein.fasta"
        
        with open(protein_fasta, 'w') as f:
            for variant in variants:
                if not self.corn_reference:
                    continue
                
                # Apply mutations
                mutated_seq = self.apply_mutations_to_sequence(
                    variant['mutations'],
                    self.corn_reference
                )
                
                # Write FASTA entry
                f.write(f">{variant['variant_id']} | {variant['description']}\n")
                f.write(f"# Mutations: {', '.join(variant['mutations'])}\n")
                f.write(f"# Priority: {variant['priority']} | Source: Level {variant['source_level']}\n")
                f.write(f"{mutated_seq}\n\n")
        
        print(f"✓ Protein sequences: {protein_fasta}")
        
        # Summary JSON
        summary_json = Path(output_dir) / "variant_summary.json"
        with open(summary_json, 'w') as f:
            json.dump(variants, f, indent=2)
        
        print(f"✓ Summary JSON: {summary_json}")
        
        # Human-readable report
        report_txt = Path(output_dir) / "variant_report.txt"
        with open(report_txt, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("CORN D1 PROTEIN ENGINEERING VARIANTS\n")
            f.write("Phylogenetically-Informed Design\n")
            f.write("=" * 80 + "\n\n")
            
            for variant in variants:
                f.write(f"\n{'='*80}\n")
                f.write(f"VARIANT: {variant['variant_id']}\n")
                f.write(f"{'='*80}\n")
                f.write(f"Description: {variant['description']}\n")
                f.write(f"Priority: {variant['priority']}\n")
                f.write(f"Source Level: {variant['source_level']}\n")
                f.write(f"Expected Improvement: {variant['expected_improvement']}\n")
                f.write(f"\nMutations ({len(variant['mutations'])}):\n")
                for mut in variant['mutations']:
                    f.write(f"  • {mut}\n")
                f.write("\n")
        
        print(f"✓ Human-readable report: {report_txt}")
        
        return {
            "protein_fasta": str(protein_fasta),
            "summary_json": str(summary_json),
            "report": str(report_txt)
        }


def main():
    """
    Main integration workflow
    """
    
    print("=" * 80)
    print("PHYLOGENETIC VARIANT GENERATION")
    print("Integrating with GardenX Pipeline")
    print("=" * 80)
    print()
    
    # Initialize generator
    generator = PhylogeneticVariantGenerator()
    
    # Load data
    if not generator.load_data():
        print("\n❌ Failed to load data. Please run phylogenetic_analysis.py first!")
        return
    
    # Rank mutations
    ranked = generator.rank_mutations()
    
    if ranked is not None:
        print(f"\n📊 Top 10 Mutations:")
        print(ranked[['position', 'reference_aa', 'mutant_aa', 'level', 'total_score']].head(10).to_string())
    
    # Generate variants
    variants = generator.generate_variant_combinations(top_n=10, max_mutations_per_variant=5)
    
    if not variants:
        print("\n❌ No variants generated")
        return
    
    # Print summary
    print(f"\n📋 VARIANT SUMMARY:")
    for v in variants:
        print(f"\n  {v['variant_id']}: {v['description']}")
        print(f"    Mutations: {', '.join(v['mutations'])}")
        print(f"    Priority: {v['priority']} | Source: Level {v['source_level']}")
    
    # Generate synthesis files
    files = generator.generate_synthesis_files(variants)
    
    print("\n" + "=" * 80)
    print("✅ VARIANT GENERATION COMPLETE")
    print("=" * 80)
    print(f"\nOutput files:")
    for file_type, filepath in files.items():
        print(f"  • {file_type}: {filepath}")
    
    print("\n🎯 NEXT STEPS:")
    print("  1. Review variant_report.txt for detailed descriptions")
    print("  2. Send corn_d1_variants_protein.fasta to synthesis company")
    print("  3. Run comparative analysis with master_pipeline_v1.py")


if __name__ == "__main__":
    main()