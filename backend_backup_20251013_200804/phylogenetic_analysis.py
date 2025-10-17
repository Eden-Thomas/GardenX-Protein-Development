# backend/phylogenetic_analysis.py
"""
Phylogenetically-informed D1 protein engineering
Starting from most closely related phyollogeny → expanding outward
"""

from Bio import Entrez, AlignIO, Phylo
from Bio.Align.Applications import ClustalwCommandline
import pandas as pd
from typing import Dict, List, Tuple

Entrez.email = "tgpaveglio@utexas.edu"

class PhylogeneticD1Analysis:
    """
    Progressive analysis from close to distant relatives
    """
    
    def __init__(self):
        self.phylogenetic_levels = {
            "level_1_intraspecific": {
                "description": "Corn varieties and wild relatives",
                "evolutionary_distance": 0.001,  # ~100k years
                "sequence_identity": 0.99,
                "organisms": [
                    ("Zea mays", "B73", "reference"),
                    ("Zea mays", "tropical_landrace", "heat_adapted"),
                    ("Zea mays ssp. parviglumis", "teosinte", "wild_ancestor"),
                ]
            },
            "level_2_close_relatives": {
                "description": "Close C4 grass relatives",
                "evolutionary_distance": 12,  # million years
                "sequence_identity": 0.90,
                "organisms": [
                    ("Sorghum bicolor", "sorghum", "heat_drought_adapted"),
                    ("Saccharum officinarum", "sugarcane", "tropical"),
                    ("Miscanthus", "miscanthus", "high_biomass"),
                ]
            },
            "level_3_distant_monocots": {
                "description": "Other grasses (C3 photosynthesis)",
                "evolutionary_distance": 50,
                "sequence_identity": 0.80,
                "organisms": [
                    ("Triticum aestivum", "wheat", "temperate"),
                    ("Oryza sativa", "rice", "aquatic"),
                    ("Hordeum vulgare", "barley", "temperate"),
                ]
            },
            "level_4_very_distant": {
                "description": "Distant monocots",
                "evolutionary_distance": 120,
                "sequence_identity": 0.70,
                "organisms": [
                    ("Phoenix dactylifera", "date_palm", "extreme_heat"),
                ]
            },
            "level_5_dicots": {
                "description": "Dicot outgroup",
                "evolutionary_distance": 170,
                "sequence_identity": 0.65,
                "organisms": [
                    ("Arabidopsis thaliana", "arabidopsis", "model_plant"),
                    ("Glycine max", "soybean", "legume"),
                ]
            },
            "level_6_extremophiles": {
                "description": "Cyanobacteria extremophiles",
                "evolutionary_distance": 2000,
                "sequence_identity": 0.50,
                "organisms": [
                    ("Thermosynechococcus vulcanus", "thermophile", "57C"),
                    ("Thermosynechococcus elongatus", "thermophile", "55C"),
                ]
            }
        }
    
    def download_level_1_corn_variants(self) -> Dict[str, str]:
        """
        Download D1 sequences from different corn varieties
        Focus on heat-adapted tropical varieties
        """
        
        print("🌽 LEVEL 1: Downloading Corn D1 Variants\n")
        
        sequences = {}
        
        # Search strategies for corn variants
        search_terms = [
            "Zea mays[Organism] AND psbA[Gene] AND (tropical OR heat)",
            "Zea mays[Organism] AND D1 protein[Protein Name]",
            "Zea mays subsp. parviglumis[Organism] AND psbA",
            "Zea diploperennis[Organism] AND psbA",
        ]
        
        for search in search_terms:
            try:
                handle = Entrez.esearch(
                    db="protein",
                    term=search,
                    retmax=5
                )
                results = Entrez.read(handle)
                handle.close()
                
                for seq_id in results['IdList']:
                    handle = Entrez.efetch(
                        db="protein",
                        id=seq_id,
                        rettype="fasta",
                        retmode="text"
                    )
                    seq_data = handle.read()
                    handle.close()
                    
                    # Parse header
                    header = seq_data.split('\n')[0]
                    sequences[seq_id] = {
                        "sequence": seq_data,
                        "header": header,
                        "level": 1
                    }
                    
                    print(f"✓ Downloaded: {header[:60]}...")
                    
            except Exception as e:
                print(f"⚠ Error with search '{search}': {e}")
        
        print(f"\n✅ Level 1 Complete: {len(sequences)} corn variants downloaded\n")
        return sequences
    
    def download_level_2_close_relatives(self) -> Dict[str, str]:
        """
        Download close C4 grass relatives (Sorghum, Sugarcane, Miscanthus)
        """
        
        print("🌾 LEVEL 2: Downloading Close C4 Grass Relatives\n")
        
        sequences = {}
        
        organisms = [
            "Sorghum bicolor",
            "Saccharum officinarum", 
            "Miscanthus"
        ]
        
        for organism in organisms:
            try:
                search = f"{organism}[Organism] AND psbA[Gene]"
                handle = Entrez.esearch(db="protein", term=search, retmax=3)
                results = Entrez.read(handle)
                handle.close()
                
                if results['IdList']:
                    seq_id = results['IdList'][0]  # Get best match
                    handle = Entrez.efetch(
                        db="protein",
                        id=seq_id,
                        rettype="fasta",
                        retmode="text"
                    )
                    seq_data = handle.read()
                    handle.close()
                    
                    sequences[organism] = {
                        "sequence": seq_data,
                        "level": 2,
                        "note": "Close C4 relative"
                    }
                    
                    print(f"✓ {organism}: Downloaded")
                
            except Exception as e:
                print(f"⚠ {organism}: {e}")
        
        print(f"\n✅ Level 2 Complete: {len(sequences)} close relatives\n")
        return sequences
    
    def analyze_progressive_mutations(self, sequences: Dict) -> pd.DataFrame:
        """
        Identify mutations progressively across phylogenetic levels
        """
        
        print("🔬 ANALYZING PROGRESSIVE MUTATIONS\n")
        
        # Organize by level
        by_level = {}
        for seq_id, data in sequences.items():
            level = data.get('level', 0)
            if level not in by_level:
                by_level[level] = []
            by_level[level].append(data['sequence'])
        
        mutations = []
        
        # For each level, compare to reference corn
        reference = None  # Will be B73 corn
        
        for level in sorted(by_level.keys()):
            if level == 1 and not reference:
                # First corn sequence is reference
                reference = by_level[level][0]
                continue
            
            print(f"\n📊 Level {level} mutations vs. corn reference:")
            
            for seq in by_level[level]:
                # Align and find differences
                differences = self.find_mutations(reference, seq)
                
                for pos, ref_aa, mut_aa in differences:
                    mutations.append({
                        "level": level,
                        "position": pos,
                        "reference_aa": ref_aa,
                        "mutant_aa": mut_aa,
                        "conservation_score": self.calculate_conservation(pos, by_level)
                    })
                    
                    print(f"  Position {pos}: {ref_aa}→{mut_aa}")
        
        return pd.DataFrame(mutations)
    
    def find_mutations(self, seq1: str, seq2: str) -> List[Tuple[int, str, str]]:
        """Simple pairwise comparison to find mutations"""
        
        # Extract just sequence (remove FASTA header)
        s1 = ''.join([line for line in seq1.split('\n')[1:]])
        s2 = ''.join([line for line in seq2.split('\n')[1:]])
        
        mutations = []
        min_len = min(len(s1), len(s2))
        
        for i in range(min_len):
            if s1[i] != s2[i] and s1[i] != '-' and s2[i] != '-':
                mutations.append((i+1, s1[i], s2[i]))
        
        return mutations
    
    def calculate_conservation(self, position: int, sequences_by_level: Dict) -> float:
        """
        Calculate how conserved a position is across levels
        Higher score = more conserved = safer to mutate
        """
        # Simplified - would use proper MSA in production
        return 0.85  # Placeholder
    
    def generate_recommendations(self, mutations_df: pd.DataFrame) -> List[Dict]:
        """
        Recommend mutations based on phylogenetic distance
        """
        
        print("\n🎯 MUTATION RECOMMENDATIONS (Ordered by Safety)\n")
        
        recommendations = []
        
        # Priority 1: Mutations found in Level 1 (corn varieties)
        level1_muts = mutations_df[mutations_df['level'] == 1]
        if not level1_muts.empty:
            print("⭐ HIGHEST PRIORITY - Found in other corn varieties:")
            for _, row in level1_muts.iterrows():
                print(f"   Position {row['position']}: {row['reference_aa']}→{row['mutant_aa']}")
                recommendations.append({
                    "position": row['position'],
                    "mutation": f"{row['reference_aa']}{row['position']}{row['mutant_aa']}",
                    "priority": "HIGHEST",
                    "source": "Natural corn variation",
                    "risk": "Very Low"
                })
        
        # Priority 2: Level 2 (close C4 relatives)
        level2_muts = mutations_df[mutations_df['level'] == 2]
        if not level2_muts.empty:
            print("\n⭐⭐ HIGH PRIORITY - Found in close C4 relatives (Sorghum, etc):")
            for _, row in level2_muts.head(5).iterrows():
                print(f"   Position {row['position']}: {row['reference_aa']}→{row['mutant_aa']}")
                recommendations.append({
                    "position": row['position'],
                    "mutation": f"{row['reference_aa']}{row['position']}{row['mutant_aa']}",
                    "priority": "HIGH",
                    "source": "Close C4 grass relatives",
                    "risk": "Low"
                })
        
        return recommendations


def main():
    """
    Run progressive phylogenetic analysis
    """
    
    print("=" * 70)
    print("PHYLOGENETICALLY-INFORMED D1 PROTEIN ENGINEERING")
    print("Progressive Analysis: Corn → Relatives → Extremophiles")
    print("=" * 70)
    print()
    
    analyzer = PhylogeneticD1Analysis()
    
    # Download progressively
    all_sequences = {}
    
    # Level 1: Corn variants (START HERE)
    level1 = analyzer.download_level_1_corn_variants()
    all_sequences.update(level1)
    
    # Level 2: Close relatives
    level2 = analyzer.download_level_2_close_relatives()
    all_sequences.update(level2)
    
    # Analyze mutations
    mutations_df = analyzer.analyze_progressive_mutations(all_sequences)
    
    # Generate recommendations
    recommendations = analyzer.generate_recommendations(mutations_df)
    
    # Save results
    mutations_df.to_csv("data/phylogenetic_mutations.csv", index=False)
    
    print("\n" + "=" * 70)
    print(f"✅ ANALYSIS COMPLETE")
    print(f"   Total sequences analyzed: {len(all_sequences)}")
    print(f"   Mutations identified: {len(mutations_df)}")
    print(f"   High-confidence recommendations: {len(recommendations)}")
    print("=" * 70)


if __name__ == "__main__":
    main()