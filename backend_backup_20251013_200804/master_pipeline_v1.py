"""
master_pipeline_v1.py
Complete V1 D1 Protein Engineering Pipeline
Integrates all components: Data, Structure, Function, Combinations, and Evo 2
"""
import os
from pathlib import Path
from dotenv import load_dotenv


# Loading environment variables from .env file
env_path = Path(__file__).parent.parent / '.env'
load_dotenv(dotenv_path=env_path)

# Tag API keys from environment
NEUROSNAP_API_KEY = os.getenv("NEUROSNAP_API_KEY")
NVIDIA_API_KEY = os.getenv("NVIDIA_API_KEY")

# Importing all components
import numpy as np
from typing import Dict, List, Optional
import json
from pathlib import Path
from datetime import datetime

from data_collection import ExperimentalDataManager, RealDataMLTrainer
from structure_aware_design import StructureGuidedMutationDesigner
from functional_design import FunctionalValidator, CombinatorialDesigner, ExperimentalPrioritizer
from evo2_integration import Evo2D1Pipeline
from neurosnap_api import HybridPredictionPipeline


class D1EngineeringPipelineV1:
    """
    Complete V1 Pipeline for D1 Protein Engineering
    
    Pipeline Steps:
    1. Data Collection - Load experimental data
    2. Comparative Analysis - Learn from thermophiles
    3. ML Prediction - Predict stability changes
    4. Structure Filtering - Protect critical sites
    5. Functional Validation - Ensure photosynthesis preserved
    6. Combinatorial Design - Generate multi-mutation variants
    7. Evo 2 Validation - DNA-level validation
    8. Prioritization - Rank for experimental testing
    """
    
    def __init__(self,
                 neurosnap_api_key: Optional[str] = None,
                 evo2_api_key: Optional[str] = None):
        """
        Initialize complete V1 pipeline
        
        Args:
            neurosnap_api_key: NeuroSnap API key (AlphaFold2, NeuroFold)
            evo2_api_key: NVIDIA API key for Evo 2
        """
        
        print("="*70)
        print("INITIALIZING D1 PROTEIN ENGINEERING PIPELINE V1")
        print("="*70)
        
        # Component 1: Data Management
        print("\n[1/7] Loading experimental data system...")
        self.data_manager = ExperimentalDataManager()
        self.data_trainer = RealDataMLTrainer(self.data_manager)
        
        # Component 2: Comparative Analysis (existing)
        print("[2/7] Loading comparative analysis engine...")
        from analysis_engine import ComparativeAnalyzer
        self.comparative_analyzer = ComparativeAnalyzer()
        
        # Component 3: Structure-Aware Design
        print("[3/7] Loading structure-aware design system...")
        self.structure_designer = StructureGuidedMutationDesigner()
        
        # Component 4: Functional Validation
        print("[4/7] Loading functional validation system...")
        self.functional_validator = FunctionalValidator()
        self.combinatorial_designer = CombinatorialDesigner()
        
        # Component 5: NeuroSnap Integration
        print("[5/7] Connecting to NeuroSnap API...")
        if neurosnap_api_key:
            self.neurosnap_pipeline = HybridPredictionPipeline(neurosnap_api_key)
            self.has_neurosnap = True
        else:
            self.has_neurosnap = False
            print("   Warning: No NeuroSnap API key - AlphaFold validation disabled")
        
        # Component 6: Evo 2 Integration
        print("[6/7] Connecting to Evo 2...")
        if evo2_api_key:
            self.evo2_pipeline = Evo2D1Pipeline(evo2_api_key)
            self.has_evo2 = True
        else:
            self.has_evo2 = False
            print("   Warning: No Evo 2 API key - DNA-level validation disabled")
        
        # Component 7: Experiment Prioritization
        print("[7/7] Loading experiment prioritization...")
        self.exp_prioritizer = ExperimentalPrioritizer()
        
        print("\n✓ Pipeline initialization complete!")
        print(f"  - NeuroSnap: {'ENABLED' if self.has_neurosnap else 'DISABLED'}")
        print(f"  - Evo 2: {'ENABLED' if self.has_evo2 else 'DISABLED'}")
    
    def run_complete_analysis(self,
                             reference_sequence: str,
                             comparison_species: List[str],
                             organism: str = 'corn',
                             output_dir: str = 'results_v1') -> Dict:
        """
        Run complete V1 analysis pipeline
        
        Args:
            reference_sequence: D1 protein sequence to optimize
            comparison_species: Species to compare (e.g., ['wheat', 'rice', 'thermo'])
            organism: Target organism for optimization
            output_dir: Directory to save results
        
        Returns:
            Complete analysis results dictionary
        """
        
        # Setup output directory
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        run_id = f"d1_analysis_{timestamp}"
        
        print("\n" + "="*70)
        print(f"RUNNING COMPLETE V1 ANALYSIS - {run_id}")
        print("="*70)
        
        results = {
            'run_id': run_id,
            'timestamp': timestamp,
            'reference_sequence': reference_sequence,
            'organism': organism,
            'pipeline_stages': {}
        }
        
        # ===================================================================
        # STAGE 1: COMPARATIVE ANALYSIS
        # ===================================================================
        print("\n" + "="*70)
        print("STAGE 1: COMPARATIVE GENOMICS ANALYSIS")
        print("="*70)
        
        from analysis_engine import run_comprehensive_analysis
        
        comparative_results = run_comprehensive_analysis(
            comparison_species,
            reference_sequence
        )
        
        results['pipeline_stages']['comparative'] = comparative_results
        
        print(f"\n✓ Stage 1 Complete:")
        print(f"  - Generated {len(comparative_results['recommendations'])} initial candidates")
        print(f"  - Conservation analysis: {comparative_results['conservation']['conserved_positions']} conserved positions")
        
        # ===================================================================
        # STAGE 2: RETRAIN ML WITH REAL DATA (if available)
        # ===================================================================
        print("\n" + "="*70)
        print("STAGE 2: ML MODEL TRAINING WITH REAL DATA")
        print("="*70)
        
        data_stats = self.data_manager.calculate_data_statistics()
        
        if data_stats.get('total_mutations', 0) > 0:
            print(f"Found {data_stats['total_mutations']} experimental data points")
            
            if data_stats['total_mutations'] < 10:
                print(f"⚠️  Warning: Less than 10 data points. Results will be unreliable.")
                print(f"   Please add more experimental data to improve predictions.")
            
            try:
                training_results = self.data_trainer.train_on_real_data()
                
                if training_results and isinstance(training_results, dict):
                    results['pipeline_stages']['ml_training'] = training_results
                    print(f"✓ Stage 2 Complete:")
                    print(f"  - Best model: {training_results.get('best_model', 'N/A')}")
                    
                    best_model = training_results.get('best_model')
                    if best_model and 'model_performance' in training_results:
                        model_perf = training_results['model_performance'].get(best_model, {})
                        r2 = model_perf.get('r2', 0)
                        print(f"  - R²: {r2:.3f}")
                    else:
                        print(f"  - R²: N/A (insufficient data)")
                else:
                    print(f"⚠️  Insufficient data for reliable model training")
                    print(f"  - Using comparative genomics predictions only")
                    results['pipeline_stages']['ml_training'] = {
                        'mode': 'insufficient_data',
                        'note': 'Less than minimum required data points'
                    }
                    
            except Exception as e:
                print(f"⚠️  Model training failed: {str(e)}")
                print(f"  - Continuing with comparative genomics only")
                results['pipeline_stages']['ml_training'] = {
                    'mode': 'failed',
                    'error': str(e)
                }
        else:
            print("No experimental data available - using comparative genomics only")
            print("  Add data to data/experimental_thermostability.csv to improve predictions")
            results['pipeline_stages']['ml_training'] = {
                'mode': 'no_data',
                'note': 'No experimental data found'
            }
        
        # ===================================================================
        # STAGE 3: STRUCTURE-AWARE FILTERING
        # ===================================================================
        print("\n" + "="*70)
        print("STAGE 3: STRUCTURE-AWARE FILTERING")
        print("="*70)
        
        if self.has_neurosnap:
            try:
                print("Predicting structure with AlphaFold2...")
                structure_result = self.neurosnap_pipeline.neurosnap.predict_structure(
                    reference_sequence,
                    job_name=f"{organism}_d1_reference"
                )
                
                pdb_content = structure_result['pdb_structure']
                structure_data = self.structure_designer.analyzer.parse_pdb_structure(pdb_content)
                
                # Filter mutations
                structure_filtered = self.structure_designer.filter_mutations_by_structure(
                    comparative_results['recommendations'],
                    structure_data
                )
                
                results['pipeline_stages']['structure_filtering'] = {
                    'n_input': len(comparative_results['recommendations']),
                    'n_output': len(structure_filtered),
                    'n_rejected': len(comparative_results['recommendations']) - len(structure_filtered),
                    'structure_confidence': structure_result['confidence']
                }
                
                print(f"✓ Stage 3 Complete:")
                print(f"  - Filtered {len(comparative_results['recommendations'])} → {len(structure_filtered)} mutations")
                print(f"  - AlphaFold confidence: {structure_result['confidence']:.2f}")
                
            except Exception as e:
                print(f"⚠️  Structure prediction failed: {str(e)}")
                print(f"  - Continuing without structure filtering")
                structure_filtered = comparative_results['recommendations']
                results['pipeline_stages']['structure_filtering'] = {
                    'mode': 'failed',
                    'error': str(e)
                }
            
        else:
            print("Structure filtering disabled (no NeuroSnap API)")
            print("  - Continuing with all comparative genomics predictions")
            structure_filtered = comparative_results['recommendations']
            results['pipeline_stages']['structure_filtering'] = {'mode': 'disabled'}
        
        # ===================================================================
        # STAGE 4: FUNCTIONAL VALIDATION
        # ===================================================================
        print("\n" + "="*70)
        print("STAGE 4: FUNCTIONAL VALIDATION")
        print("="*70)
        
        functional_mutations = self.functional_validator.filter_by_function(
            structure_filtered,
            min_function_score=0.7
        )
        
        results['pipeline_stages']['functional'] = {
            'n_input': len(structure_filtered),
            'n_output': len(functional_mutations),
            'n_rejected_functional': len(structure_filtered) - len(functional_mutations)
        }
        
        print(f"✓ Stage 4 Complete:")
        print(f"  - Kept {len(functional_mutations)}/{len(structure_filtered)} functional mutations")
        
        # ===================================================================
        # STAGE 5: COMBINATORIAL DESIGN
        # ===================================================================
        print("\n" + "="*70)
        print("STAGE 5: COMBINATORIAL VARIANT GENERATION")
        print("="*70)
        
        all_variants = self.combinatorial_designer.generate_combinations(
            functional_mutations,
            max_mutations=3,
            min_spacing=15
        )
        
        ranked_variants = self.combinatorial_designer.rank_variants(
            all_variants,
            reference_sequence
        )
        
        results['pipeline_stages']['combinatorial'] = {
            'n_singles': len([v for v in ranked_variants if v['n_mutations'] == 1]),
            'n_doubles': len([v for v in ranked_variants if v['n_mutations'] == 2]),
            'n_triples': len([v for v in ranked_variants if v['n_mutations'] == 3]),
            'total_variants': len(ranked_variants),
            'top_variant': ranked_variants[0] if ranked_variants else None
        }
        
        print(f"✓ Stage 5 Complete:")
        print(f"  - Generated {len(ranked_variants)} total variants")
        print(f"  - Singles: {results['pipeline_stages']['combinatorial']['n_singles']}")
        print(f"  - Doubles: {results['pipeline_stages']['combinatorial']['n_doubles']}")
        print(f"  - Triples: {results['pipeline_stages']['combinatorial']['n_triples']}")
        
        # ===================================================================
        # STAGE 6: EVO 2 DNA-LEVEL VALIDATION
        # ===================================================================
        print("\n" + "="*70)
        print("STAGE 6: EVO 2 DNA-LEVEL VALIDATION")
        print("="*70)
        
        if self.has_evo2:
            try:
                print("Running Evo 2 validation on top candidates...")
                
                # Take top 20 variants for Evo 2 validation
                top_variants = ranked_variants[:20]
                
                evo2_results = self.evo2_pipeline.engineer_d1_gene(
                    reference_sequence,
                    functional_mutations[:20],  # Top 20 single mutations
                    organism=organism
                )
                
                results['pipeline_stages']['evo2'] = {
                    'n_tested': len(top_variants),
                    'n_approved': evo2_results['n_mutations_applied'],
                    'approval_rate': evo2_results['evo2_validation']['approval_rate'],
                    'gene_construct_ready': evo2_results['ready_for_lab']
                }
                
                print(f"✓ Stage 6 Complete:")
                print(f"  - Approved {evo2_results['n_mutations_applied']}/{len(top_variants)} variants")
                print(f"  - Gene construct ready: {evo2_results['ready_for_lab']}")
                
                # Update variants with Evo 2 scores
                evo2_approved_positions = set(m['position'] for m in evo2_results['mutations'])
                
                for variant in ranked_variants:
                    variant_positions = set(m['position'] for m in variant['mutations'])
                    if variant_positions.issubset(evo2_approved_positions):
                        variant['evo2_validated'] = True
                    else:
                        variant['evo2_validated'] = False
                        
            except Exception as e:
                print(f"⚠️  Evo 2 validation failed: {str(e)}")
                print(f"  - Continuing without DNA-level validation")
                results['pipeline_stages']['evo2'] = {
                    'mode': 'failed',
                    'error': str(e)
                }
            
        else:
            print("Evo 2 validation disabled (no API key)")
            print("  - Skipping DNA-level validation")
            results['pipeline_stages']['evo2'] = {'mode': 'disabled'}
        
        # ===================================================================
        # STAGE 7: EXPERIMENTAL PRIORITIZATION
        # ===================================================================
        print("\n" + "="*70)
        print("STAGE 7: EXPERIMENTAL TESTING PRIORITIZATION")
        print("="*70)
        
        experimental_plan = self.exp_prioritizer.create_experimental_plan(
            ranked_variants,
            budget=20
        )
        
        # Export sequences
        synthesis_file = self.exp_prioritizer.export_sequences_for_synthesis(
            experimental_plan,
            output_file=str(output_path / f"{run_id}_synthesis.fasta")
        )
        
        results['pipeline_stages']['prioritization'] = {
            'wave_1_count': len(experimental_plan['wave_1_singles']),
            'wave_2_count': len(experimental_plan['wave_2_doubles']),
            'wave_3_count': len(experimental_plan['wave_3_triples']),
            'synthesis_file': synthesis_file
        }
        
        print(f"✓ Stage 7 Complete:")
        print(f"  - Wave 1 (singles): {len(experimental_plan['wave_1_singles'])} variants")
        print(f"  - Wave 2 (doubles): {len(experimental_plan['wave_2_doubles'])} variants")
        print(f"  - Wave 3 (triples): {len(experimental_plan['wave_3_triples'])} variants")
        print(f"  - Synthesis file: {synthesis_file}")
        
        # ===================================================================
        # FINAL RESULTS & EXPORT
        # ===================================================================
        print("\n" + "="*70)
        print("GENERATING FINAL RESULTS")
        print("="*70)
        
        results['final_recommendations'] = {
            'all_variants': ranked_variants,
            'experimental_plan': experimental_plan,
            'top_10_variants': ranked_variants[:10],
            'ready_for_synthesis': [
                v for v in ranked_variants 
                if v.get('evo2_validated', True) and v['overall_score'] > 0.7
            ]
        }
        
        # Save complete results
        results_file = output_path / f"{run_id}_complete_results.json"
        with open(results_file, 'w') as f:
            # Convert numpy types to native Python for JSON
            json_results = self._prepare_for_json(results)
            json.dump(json_results, f, indent=2)
        
        print(f"\n✓ Results saved to: {results_file}")
        
        # Generate summary report
        self._generate_summary_report(results, output_path / f"{run_id}_summary.txt")
        
        print("\n" + "="*70)
        print("V1 ANALYSIS COMPLETE!")
        print("="*70)
        print(f"\nRun ID: {run_id}")
        print(f"Results directory: {output_path}")
        print(f"\nFiles generated:")
        print(f"  - {run_id}_complete_results.json (full data)")
        print(f"  - {run_id}_summary.txt (human-readable report)")
        print(f"  - {run_id}_synthesis.fasta (sequences for gene synthesis)")
        
        return results
    
    def _prepare_for_json(self, obj):
        """Convert numpy types to native Python for JSON serialization"""
        # Handle numpy numeric types
        if isinstance(obj, (np.integer, np.int64, np.int32, np.int16, np.int8)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64, np.float32, np.float16)):
            return float(obj)
        # Handle numpy boolean types - THIS IS THE KEY FIX
        elif isinstance(obj, np.bool_):
            return bool(obj)
        # Handle numpy arrays
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        # Handle dictionaries recursively
        elif isinstance(obj, dict):
            return {key: self._prepare_for_json(value) for key, value in obj.items()}
        # Handle lists recursively
        elif isinstance(obj, list):
            return [self._prepare_for_json(item) for item in obj]
        # Handle tuples (convert to list)
        elif isinstance(obj, tuple):
            return [self._prepare_for_json(item) for item in obj]
        # Handle None and basic Python types
        elif obj is None or isinstance(obj, (str, int, float, bool)):
            return obj
        # Fallback for unknown types
        else:
            try:
                return str(obj)
            except:
                return None
    
    def _generate_summary_report(self, results: Dict, output_file: Path):
        """Generate human-readable summary report"""
        
        with open(output_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("D1 PROTEIN ENGINEERING - V1 ANALYSIS SUMMARY\n")
            f.write("="*70 + "\n\n")
            
            f.write(f"Run ID: {results['run_id']}\n")
            f.write(f"Timestamp: {results['timestamp']}\n")
            f.write(f"Organism: {results['organism']}\n\n")
            
            f.write("PIPELINE STAGES SUMMARY\n")
            f.write("-"*70 + "\n\n")
            
            # Stage summaries
            for stage_name, stage_data in results['pipeline_stages'].items():
                f.write(f"{stage_name.upper()}:\n")
                if isinstance(stage_data, dict):
                    for key, value in stage_data.items():
                        f.write(f"  {key}: {value}\n")
                f.write("\n")
            
            # Top recommendations
            f.write("TOP 10 RECOMMENDED VARIANTS\n")
            f.write("-"*70 + "\n\n")
            
            for i, variant in enumerate(results['final_recommendations']['top_10_variants'], 1):
                f.write(f"{i}. {variant['variant_id']}\n")
                f.write(f"   Predicted ΔΔG: {variant['predicted_ddg']:.2f} kcal/mol\n")
                f.write(f"   Estimated ΔTm: {variant['estimated_tm_increase']:.1f} °C\n")
                f.write(f"   Functional Score: {variant['functional_score']:.2f}\n")
                f.write(f"   Overall Score: {variant['overall_score']:.3f}\n")
                f.write(f"   Recommendation: {variant['recommendation']}\n")
                f.write("\n")
        
        print(f"✓ Summary report saved to: {output_file}")


# Example usage script
def main():
    """Example usage of complete V1 pipeline"""
    
    # API keys (get from providers)
    NEUROSNAP_KEY = "your_neurosnap_api_key"  # From neurosnap.ai
    EVO2_KEY = "your_nvidia_api_key"  # From NVIDIA
    
    # Initialize pipeline
    pipeline = D1EngineeringPipelineV1(
        neurosnap_api_key=NEUROSNAP_KEY,
        evo2_api_key=EVO2_KEY
    )
    
    # Reference D1 sequence (corn)
    reference_d1 = """
    MTAILERRESESLWGRFCNWITGTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDID
    GIREPVSGSLLYGNNIITGAVVPSSNAILVHFYPIWEAASLDEWLYNGGPYQLIIFHFLLG
    ASCYMGREWELSFRLGMRPWICVAYSAPLASAFAVFLIYPIGQGSFSDGMPLGISGTFNFM
    IVFQAEHNILMHPFHMLGVAGVFGGSLLFSAMHGSLVTSSLIRETTENESANEGYRFGQEE
    ETYNIVAHAGYGRLLFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQS
    VVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAAEG
    """.replace('\n', '').replace(' ', '')
    
    # Comparison species
    comparison_species = ['wheat', 'rice', 'thermo']
    
    # Run complete analysis
    results = pipeline.run_complete_analysis(
        reference_sequence=reference_d1,
        comparison_species=comparison_species,
        organism='corn',
        output_dir='results_v1'
    )
    
    print("\n🎉 Analysis complete! Check results_v1/ directory for outputs.")
    
    return results


if __name__ == "__main__":
    main()