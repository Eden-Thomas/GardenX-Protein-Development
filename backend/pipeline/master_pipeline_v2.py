"""
master_pipeline_v2.py
Master Pipeline - Orchestrates all 3 tracks + consensus scoring

Complete workflow:
1. Track 1: Evolutionary-Guided Design
2. Track 2: Generative AI Design
3. Track 3: Structure-Guided Design
4. Consensus Scoring & Integration
5. Final Report Generation
"""

import sys
from pathlib import Path
from datetime import datetime
import json

# Add parent directory to path
sys.path.append(str(Path(__file__).parent.parent))

from analysis_engine import REFERENCE_SEQUENCES, GROWTH_TEMPERATURES
from pipeline.track1_evolutionary import EvolutionaryGuidedDesigner
from pipeline.track2_generative import GenerativeAIDesigner
from pipeline.track3_structural import StructuralDesigner
from pipeline.consensus_scorer import ConsensusScorer


class MasterPipeline:
    """
    Master pipeline orchestrating complete variant design workflow
    """
    
    def __init__(self, api_key: str = None):
        """
        Args:
            api_key: Neurosnap API key (optional, reads from env)
        """
        self.api_key = api_key
        self.maize_sequence = REFERENCE_SEQUENCES['maize']
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.results_dir = Path(f"../data/results/master_run_{self.timestamp}")
        
        # Results storage
        self.track1_results = None
        self.track2_results = None
        self.track3_results = None
        self.consensus_results = None
    
    def run_complete_pipeline(self, 
                            # Track 1 parameters
                            track1_min_correlation: float = 0.3,
                            track1_max_conservation: float = 0.7,
                            track1_p_threshold: float = 0.10,
                            track1_n_variants: int = 20,
                            # Track 2 parameters
                            track2_n_efficient_evo: int = 15,
                            track2_n_rfdiffusion: int = 10,
                            track2_max_mutations: int = 3,
                            # Track 3 parameters
                            track3_n_variants: int = 75) -> dict:
        """
        Run complete 3-track pipeline with consensus scoring
        
        Returns:
            Dict with all results
        """
        
        print("\n" + "="*70)
        print("🌽 MAIZE D1 MASTER PIPELINE v2.0")
        print("="*70)
        print(f"Start Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Template: Zea mays D1 ({len(self.maize_sequence)} aa)")
        print(f"Results Directory: {self.results_dir}")
        print("="*70 + "\n")
        
        # Create results directory
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # ===================================================================
        # TRACK 1: EVOLUTIONARY-GUIDED DESIGN
        # ===================================================================
        print("\n" + "▶"*35)
        print("TRACK 1: EVOLUTIONARY-GUIDED DESIGN")
        print("▶"*35)
        
        track1_designer = EvolutionaryGuidedDesigner(
            sequences=REFERENCE_SEQUENCES,
            temperatures=GROWTH_TEMPERATURES
        )
        
        self.track1_results = track1_designer.run_complete_analysis(
            min_correlation=track1_min_correlation,
            max_conservation=track1_max_conservation,
            p_value_threshold=track1_p_threshold,
            n_variants=track1_n_variants
        )
        
        track1_dir = self.results_dir / "track1_evolutionary"
        track1_designer.export_results(track1_dir, self.track1_results)
        
        print(f"\n✅ Track 1 Complete: {len(self.track1_results['variants'])} variants generated")
        
        # ===================================================================
        # TRACK 2: GENERATIVE AI DESIGN
        # ===================================================================
        print("\n" + "▶"*35)
        print("TRACK 2: GENERATIVE AI DESIGN")
        print("▶"*35)
        
        track2_designer = GenerativeAIDesigner(
            maize_sequence=self.maize_sequence,
            api_key=self.api_key
        )
        
        self.track2_results = track2_designer.run_complete_analysis(
            n_efficient_evolution=track2_n_efficient_evo,
            n_rfdiffusion=track2_n_rfdiffusion,
            max_mutations=track2_max_mutations
        )
        
        track2_dir = self.results_dir / "track2_generative"
        track2_designer.export_results(track2_dir, self.track2_results)
        
        print(f"\n✅ Track 2 Complete: {len(self.track2_results['variants'])} variants generated")
        
        # ===================================================================
        # TRACK 3: STRUCTURE-GUIDED DESIGN
        # ===================================================================
        print("\n" + "▶"*35)
        print("TRACK 3: STRUCTURE-GUIDED DESIGN")
        print("▶"*35)
        
        track3_designer = StructuralDesigner(
            maize_sequence=self.maize_sequence,
            api_key=self.api_key
        )
        
        self.track3_results = track3_designer.run_complete_analysis(
            n_variants=track3_n_variants
        )
        
        track3_dir = self.results_dir / "track3_structural"
        track3_designer.export_results(track3_dir, self.track3_results)
        
        print(f"\n✅ Track 3 Complete: {len(self.track3_results['variants'])} variants generated")
        
        # ===================================================================
        # CONSENSUS SCORING & INTEGRATION
        # ===================================================================
        print("\n" + "▶"*35)
        print("CONSENSUS SCORING & INTEGRATION")
        print("▶"*35)
        
        scorer = ConsensusScorer(maize_sequence=self.maize_sequence)
        
        self.consensus_results = scorer.score_all_variants(
            self.track1_results,
            self.track2_results,
            self.track3_results
        )
        
        consensus_dir = self.results_dir / "consensus"
        scorer.export_results(consensus_dir, self.consensus_results)
        
        print(f"\n✅ Consensus Scoring Complete")
        
        # ===================================================================
        # GENERATE MASTER REPORT
        # ===================================================================
        print("\n" + "▶"*35)
        print("GENERATING MASTER REPORT")
        print("▶"*35)
        
        master_summary = self._generate_master_summary()
        
        # Export master summary
        summary_file = self.results_dir / "MASTER_SUMMARY.json"
        with open(summary_file, 'w') as f:
            json.dump(master_summary, f, indent=2, default=str)
        
        # Generate README
        self._generate_readme(master_summary)
        
        print(f"\n💾 Master summary: {summary_file}")
        print(f"💾 README: {self.results_dir / 'README.md'}")
        
        # ===================================================================
        # FINAL SUMMARY
        # ===================================================================
        print("\n" + "="*70)
        print("🎉 MASTER PIPELINE COMPLETE!")
        print("="*70)
        print(f"⏱️  Total Runtime: {master_summary['runtime']}")
        print(f"\n📊 RESULTS SUMMARY:")
        print(f"   Total Variants Generated: {master_summary['total_variants']}")
        print(f"   - Track 1 (Evolutionary):  {master_summary['track1_variants']}")
        print(f"   - Track 2 (Generative):    {master_summary['track2_variants']}")
        print(f"   - Track 3 (Structural):    {master_summary['track3_variants']}")
        print(f"\n🏆 CONSENSUS:")
        print(f"   Consensus Mutations: {master_summary['consensus_mutations']}")
        print(f"   Top Variant: {master_summary['top_variant']}")
        print(f"   Top Score: {master_summary['top_score']:.2f}")
        print(f"\n📁 All results saved to: {self.results_dir}")
        print("="*70 + "\n")
        
        return master_summary
    
    def _generate_master_summary(self) -> dict:
        """Generate comprehensive summary of all results"""
        
        end_time = datetime.now()
        start_time_str = self.timestamp
        start_time = datetime.strptime(start_time_str, "%Y%m%d_%H%M%S")
        runtime = str(end_time - start_time)
        
        summary = {
            'pipeline_version': '2.0',
            'run_timestamp': self.timestamp,
            'runtime': runtime,
            
            # Variant counts
            'total_variants': self.consensus_results['total_variants'],
            'track1_variants': len(self.track1_results['variants']),
            'track2_variants': len(self.track2_results['variants']),
            'track3_variants': len(self.track3_results['variants']),
            
            # Consensus results
            'consensus_mutations': len(self.consensus_results['consensus_mutations']),
            'consensus_mutation_list': list(self.consensus_results['consensus_mutations'].keys()),
            
            # Top variant
            'top_variant': (self.consensus_results['top_variants'][0]['variant_id'] 
                          if self.consensus_results['top_variants'] else 'N/A'),
            'top_score': (self.consensus_results['top_variants'][0]['final_score']
                        if self.consensus_results['top_variants'] else 0),
            'top_mutations': (self.consensus_results['top_variants'][0]['mutations']
                            if self.consensus_results['top_variants'] else []),
            
            # Statistics
            'score_statistics': self.consensus_results['statistics'],
            
            # Top 10 variants
            'top_10_variants': [
                {
                    'variant_id': v['variant_id'],
                    'mutations': v['mutations'],
                    'score': v['final_score'],
                    'track': v['source_track'],
                    'has_consensus': v['has_consensus']
                }
                for v in self.consensus_results['top_variants'][:10]
            ],
            
            # Directories
            'results_directory': str(self.results_dir),
            'track1_dir': str(self.results_dir / "track1_evolutionary"),
            'track2_dir': str(self.results_dir / "track2_generative"),
            'track3_dir': str(self.results_dir / "track3_structural"),
            'consensus_dir': str(self.results_dir / "consensus"),
        }
        
        return summary
    
    def _generate_readme(self, summary: dict):
        """Generate human-readable README"""
        
        readme_content = f"""# Maize D1 Variant Design - Pipeline Run Results

**Run Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}  
**Pipeline Version:** {summary['pipeline_version']}  
**Runtime:** {summary['runtime']}

---

## 📊 Summary

- **Total Variants Generated:** {summary['total_variants']}
- **Consensus Mutations Found:** {summary['consensus_mutations']}
- **Top Variant Score:** {summary['top_score']:.2f}

---

## 🎯 Results by Track

### Track 1: Evolutionary-Guided Design
- **Variants Generated:** {summary['track1_variants']}
- **Method:** Comparative analysis of heat-adapted plants and thermophiles
- **Results:** `track1_evolutionary/`

### Track 2: Generative AI Design
- **Variants Generated:** {summary['track2_variants']}
- **Methods:** Efficient Evolution + RFdiffusion
- **Results:** `track2_generative/`

### Track 3: Structure-Guided Design
- **Variants Generated:** {summary['track3_variants']}
- **Methods:** Boltz-1 structure prediction + ProteinMPNN
- **Results:** `track3_structural/`

---

## 🏆 Top 10 Candidates

"""
        
        for i, variant in enumerate(summary['top_10_variants'], 1):
            consensus_mark = "✅" if variant['has_consensus'] else ""
            mutations_str = ", ".join(variant['mutations'])
            readme_content += f"{i}. **{variant['variant_id']}** {consensus_mark}\n"
            readme_content += f"   - Mutations: {mutations_str}\n"
            readme_content += f"   - Score: {variant['score']:.2f}\n"
            readme_content += f"   - Source: {variant['track']}\n\n"
        
        readme_content += """---

## 📁 File Structure
```
master_run_TIMESTAMP/
├── MASTER_SUMMARY.json          # Complete run summary
├── README.md                     # This file
├── track1_evolutionary/          # Track 1 results
│   ├── track1_evolutionary_variants.fasta
│   ├── track1_candidate_positions.csv
│   └── track1_analysis_results.json
├── track2_generative/            # Track 2 results
│   ├── track2_generative_variants.fasta
│   └── track2_analysis_results.json
├── track3_structural/            # Track 3 results
│   ├── track3_structural_variants.fasta
│   ├── structure_info.json
│   └── track3_analysis_results.json
└── consensus/                    # Final integrated results
    ├── final_top_variants.fasta
    ├── variant_scores.csv
    ├── consensus_mutations.json
    └── consensus_results.json
```

---

## 🔬 Next Steps

1. **Review Consensus Mutations:** Check `consensus/consensus_mutations.json`
2. **Select Candidates:** Review `consensus/variant_scores.csv`
3. **Prepare for Synthesis:** Use `consensus/final_top_variants.fasta`
4. **Wet-Lab Validation:** Test top 10-15 variants for thermostability

---

## 📊 Consensus Mutations

"""
        
        if summary['consensus_mutation_list']:
            for mutation in summary['consensus_mutation_list']:
                tracks = self.consensus_results['consensus_mutations'][mutation]
                readme_content += f"- **{mutation}** - Found in: {', '.join(tracks)}\n"
        else:
            readme_content += "*No cross-track consensus mutations detected in this run.*\n"
        
        readme_content += """
---

*Generated by GardenX Maize D1 Master Pipeline v2.0*
"""
        
        readme_file = self.results_dir / "README.md"
        with open(readme_file, 'w') as f:
            f.write(readme_content)


# Command-line interface
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='GardenX Maize D1 Master Pipeline - Complete 3-Track Variant Design'
    )
    
    parser.add_argument('--api-key', type=str, help='Neurosnap API key (optional, uses env var)')
    
    # Track 1 parameters
    parser.add_argument('--t1-correlation', type=float, default=0.3,
                       help='Track 1: Min correlation (default: 0.3)')
    parser.add_argument('--t1-conservation', type=float, default=0.7,
                       help='Track 1: Max conservation (default: 0.7)')
    parser.add_argument('--t1-pvalue', type=float, default=0.10,
                       help='Track 1: P-value threshold (default: 0.10)')
    parser.add_argument('--t1-variants', type=int, default=20,
                       help='Track 1: Number of variants (default: 20)')
    
    # Track 2 parameters
    parser.add_argument('--t2-efficient-evo', type=int, default=15,
                       help='Track 2: Efficient Evolution variants (default: 15)')
    parser.add_argument('--t2-rfdiffusion', type=int, default=10,
                       help='Track 2: RFdiffusion variants (default: 10)')
    parser.add_argument('--t2-max-mutations', type=int, default=3,
                       help='Track 2: Max mutations per variant (default: 3)')
    
    # Track 3 parameters
    parser.add_argument('--t3-variants', type=int, default=75,
                       help='Track 3: Number of variants (default: 75)')
    
    args = parser.parse_args()
    
    # Run pipeline
    pipeline = MasterPipeline(api_key=args.api_key)
    
    results = pipeline.run_complete_pipeline(
        # Track 1
        track1_min_correlation=args.t1_correlation,
        track1_max_conservation=args.t1_conservation,
        track1_p_threshold=args.t1_pvalue,
        track1_n_variants=args.t1_variants,
        # Track 2
        track2_n_efficient_evo=args.t2_efficient_evo,
        track2_n_rfdiffusion=args.t2_rfdiffusion,
        track2_max_mutations=args.t2_max_mutations,
        # Track 3
        track3_n_variants=args.t3_variants
    )