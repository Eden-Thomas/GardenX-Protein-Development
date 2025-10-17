"""
backend/pipeline/master_pipeline_v3.py
Master Pipeline V3 - Production-ready orchestrator
Coordinates all tracks, validation, and consensus scoring
"""

import asyncio
import os
import sys
from typing import List, Dict, Optional, Callable
from datetime import datetime
import json
import pandas as pd
import numpy as np

# Add parent directory to path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config import Config
from neurosnap_client import NeurosnapClient
from validation.comprehensive_validator import ComprehensiveValidator

# Import tracks
from pipeline.track1_evolutionary import Track1Evolutionary
from pipeline.track2_generative import Track2Generative
from pipeline.track3_structural import Track3Structural
from pipeline.consensus_scorer import ConsensusScorer

class MasterPipelineV3:
    """
    Production-ready master pipeline
    - Real API calls only (no mock data)
    - Complete validation pipeline
    - Cost-aware processing
    """
    
    def __init__(self, config: Config = None):
        """Initialize pipeline with configuration"""
        
        self.config = config or Config()
        
        # Initialize Neurosnap client
        self.client = NeurosnapClient(self.config.NEUROSNAP_API_KEY)
        
        # Initialize tracks
        self.track1 = Track1Evolutionary(self.config, self.client)
        self.track2 = Track2Generative(self.config, self.client)
        self.track3 = Track3Structural(self.config, self.client)
        
        # Initialize consensus scorer
        self.consensus = ConsensusScorer(self.config, self.client)
        
        # Validator initialized per run
        self.validator = None
        
        # Cost tracking
        self.total_api_calls = 0
        self.estimated_cost = 0.0
        
        # Results directory
        self.results_dir = self.config.RESULTS_DIR
        os.makedirs(self.results_dir, exist_ok=True)
    
    async def run(self,
                  wt_sequence: str,
                  target_temp: int = None,
                  n_variants: int = None,
                  validate_all: bool = False,
                  progress_callback: Optional[Callable] = None) -> Dict:
        """
        Run complete D1 engineering pipeline
        
        Args:
            wt_sequence: Wild-type D1 sequence
            target_temp: Target temperature (default from config)
            n_variants: Number of variants to generate (default from config)
            validate_all: Whether to validate all variants (expensive) or just top ones
            progress_callback: Optional callback(progress, message)
        
        Returns:
            Dictionary with variants and metadata
        """
        
        # Use defaults from config if not provided
        target_temp = target_temp or self.config.TARGET_TEMPERATURE
        n_variants = n_variants or self.config.INITIAL_VARIANT_COUNT
        
        print(f"""
╔════════════════════════════════════════════════════╗
║         D1 ENGINEERING PIPELINE V3.0                ║
║              PRODUCTION RUN                         ║
║                                                     ║
║  Target Temperature: {target_temp}°C                       ║
║  Variants to Generate: {n_variants}                       ║
║  Full Validation: {validate_all}                         ║
║  API Key: {'✓' if self.config.NEUROSNAP_API_KEY else '✗'}                         ║
╚════════════════════════════════════════════════════╝
        """)
        
        # Initialize tracking
        pipeline_start = datetime.now()
        results = {
            'start_time': pipeline_start.isoformat(),
            'config': {
                'target_temp': target_temp,
                'n_variants': n_variants,
                'wt_length': len(wt_sequence)
            },
            'stages': {},
            'api_usage': {}
        }
        
        try:
            # Initialize validator with WT sequence
            print("\n📋 INITIALIZATION")
            print("─" * 50)
            
            if progress_callback:
                await progress_callback(0.05, "Initializing validator...")
            
            self.validator = ComprehensiveValidator(self.client, wt_sequence)
            await self.validator.initialize()
            
            print(f"  ✓ Validator initialized")
            print(f"  ✓ WT Tm: {self.validator.wt_tm:.1f}°C")
            
            # PHASE 1: Variant Generation (30% progress)
            print("\n🧬 PHASE 1: VARIANT GENERATION")
            print("─" * 50)
            
            if progress_callback:
                await progress_callback(0.1, "Generating variants...")
            
            all_variants = await self._generate_all_variants(
                wt_sequence,
                target_temp,
                n_variants
            )
            
            # Combine all variants
            combined_variants = []
            for track_name, track_variants in all_variants.items():
                combined_variants.extend(track_variants)
            
            results['stages']['generation'] = {
                'total_generated': len(combined_variants),
                'by_track': {
                    'track1': len(all_variants.get('track1', [])),
                    'track2': len(all_variants.get('track2', [])),
                    'track3': len(all_variants.get('track3', []))
                }
            }
            
            print(f"\n  Total variants generated: {len(combined_variants)}")
            
            if progress_callback:
                await progress_callback(0.3, f"Generated {len(combined_variants)} variants")
            
            # PHASE 2: Initial Screening with TemStaPro (50% progress)
            print("\n🔍 PHASE 2: THERMOSTABILITY SCREENING")
            print("─" * 50)
            
            if progress_callback:
                await progress_callback(0.35, "Screening variants...")
            
            screened_variants = await self._screen_variants(combined_variants, target_temp)
            
            results['stages']['screening'] = {
                'passed': len(screened_variants),
                'filtered': len(combined_variants) - len(screened_variants)
            }
            
            print(f"  ✓ {len(screened_variants)} variants passed screening")
            
            if progress_callback:
                await progress_callback(0.5, f"Screened to {len(screened_variants)} variants")
            
            # PHASE 3: Comprehensive Validation (75% progress)
            print("\n✅ PHASE 3: COMPREHENSIVE VALIDATION")
            print("─" * 50)
            
            if progress_callback:
                await progress_callback(0.55, "Validating top variants...")
            
            # Decide how many to validate based on cost considerations
            if validate_all:
                variants_to_validate = screened_variants
            else:
                # Validate top variants only (cost-conscious)
                max_validate = min(50, len(screened_variants))
                variants_to_validate = screened_variants[:max_validate]
            
            print(f"  Validating {len(variants_to_validate)} variants...")
            
            validated_variants = await self.validator.batch_validate(
                variants_to_validate,
                max_concurrent=5  # Control API load
            )
            
            # Categorize by risk
            risk_categories = self._categorize_by_risk(validated_variants)
            
            results['stages']['validation'] = {
                'total_validated': len(validated_variants),
                'low_risk': len(risk_categories['low']),
                'medium_risk': len(risk_categories['medium']),
                'high_risk': len(risk_categories['high'])
            }
            
            print(f"  ✓ Validation complete:")
            print(f"    - Low risk: {len(risk_categories['low'])}")
            print(f"    - Medium risk: {len(risk_categories['medium'])}")
            print(f"    - High risk: {len(risk_categories['high'])}")
            
            if progress_callback:
                await progress_callback(0.75, "Applying consensus scoring...")
            
            # PHASE 4: Consensus Scoring (85% progress)
            print("\n🎯 PHASE 4: CONSENSUS SCORING")
            print("─" * 50)
            
            # Score low and medium risk variants
            scorable_variants = risk_categories['low'] + risk_categories['medium']
            
            if scorable_variants:
                scored_variants = await self.consensus.score_variants(scorable_variants)
                
                results['stages']['consensus'] = {
                    'scored': len(scored_variants),
                    'avg_score': np.mean([v.get('consensus_score', 0) for v in scored_variants])
                }
                
                print(f"  ✓ Scored {len(scored_variants)} variants")
            else:
                scored_variants = []
                print("  ⚠ No variants passed validation for scoring")
            
            if progress_callback:
                await progress_callback(0.85, "Selecting final variants...")
            
            # PHASE 5: Final Selection (95% progress)
            print("\n🏆 PHASE 5: FINAL SELECTION")
            print("─" * 50)
            
            final_variants = self._select_final_variants(scored_variants)
            
            results['stages']['selection'] = {
                'final_count': len(final_variants)
            }
            
            if progress_callback:
                await progress_callback(0.95, "Generating reports...")
            
            # PHASE 6: Output Generation (100% progress)
            print("\n📊 PHASE 6: GENERATING OUTPUTS")
            print("─" * 50)
            
            output_files = self._generate_outputs(final_variants, results)
            
            results['output_files'] = output_files
            results['final_variants'] = final_variants
            
            # Calculate runtime and costs
            pipeline_end = datetime.now()
            runtime = (pipeline_end - pipeline_start).total_seconds()
            
            results['runtime_seconds'] = runtime
            results['api_usage'] = {
                'total_calls': self.client.api_calls,
                'estimated_cost': self.client.api_costs
            }
            
            if progress_callback:
                await progress_callback(1.0, "Pipeline complete!")
            
            # Print summary
            self._print_summary(final_variants, results)
            
            return results
            
        except Exception as e:
            print(f"\n❌ Pipeline failed: {str(e)}")
            import traceback
            traceback.print_exc()
            
            results['error'] = str(e)
            results['status'] = 'failed'
            
            return results
    
    async def _generate_all_variants(self, 
                                    wt_sequence: str,
                                    target_temp: int,
                                    n_variants: int) -> Dict:
        """Generate variants from all tracks"""
        
        # Distribute variants across tracks
        # Track 2 (generative with NeuroFold) gets more allocation
        track1_count = n_variants // 3
        track2_count = (n_variants // 2) + 1  # Slightly more for NeuroFold
        track3_count = n_variants - track1_count - track2_count
        
        print(f"  Variant distribution:")
        print(f"    Track 1 (Evolutionary): {track1_count}")
        print(f"    Track 2 (Generative/NeuroFold): {track2_count}")
        print(f"    Track 3 (Structural): {track3_count}")
        
        # Run tracks in parallel
        tasks = []
        
        # Only add tasks if count > 0
        if track1_count > 0:
            tasks.append(self.track1.generate_variants(wt_sequence, track1_count))
        if track2_count > 0:
            tasks.append(self.track2.generate_variants(wt_sequence, target_temp, track2_count))
        if track3_count > 0:
            tasks.append(self.track3.generate_variants(wt_sequence, track3_count))
        
        # Execute parallel generation
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # Process results
        all_variants = {}
        track_names = []
        
        if track1_count > 0:
            track_names.append('track1')
        if track2_count > 0:
            track_names.append('track2')
        if track3_count > 0:
            track_names.append('track3')
        
        for track_name, result in zip(track_names, results):
            if isinstance(result, Exception):
                print(f"  ❌ {track_name} failed: {str(result)}")
                all_variants[track_name] = []
            else:
                all_variants[track_name] = result
                print(f"  ✓ {track_name}: {len(result)} variants generated")
        
        return all_variants
    
    async def _screen_variants(self, variants: List[Dict], target_temp: int) -> List[Dict]:
        """Screen variants using TemStaPro"""
        
        print(f"  Screening {len(variants)} variants...")
        
        # Extract unique sequences (avoid duplicates)
        unique_sequences = {}
        for v in variants:
            seq = v['sequence']
            if seq not in unique_sequences:
                unique_sequences[seq] = v
        
        print(f"    {len(unique_sequences)} unique sequences")
        
        # Run TemStaPro predictions
        sequences_list = list(unique_sequences.keys())
        thermo_results = await self.client.run_temstapro_batch(sequences_list)
        
        # Add predictions to variants
        screened = []
        for seq, thermo in zip(sequences_list, thermo_results):
            variant = unique_sequences[seq]
            variant['predicted_tm'] = thermo.get('predicted_tm', 0)
            variant['thermophilic_score'] = thermo.get('thermophilic_score', 0)
            
            # Filter by predicted Tm
            if variant['predicted_tm'] >= target_temp - 5:
                screened.append(variant)
        
        # Sort by predicted Tm
        screened.sort(key=lambda x: x.get('predicted_tm', 0), reverse=True)
        
        print(f"    Best predicted Tm: {screened[0]['predicted_tm']:.1f}°C" if screened else "    No variants passed")
        
        return screened
    
    def _categorize_by_risk(self, variants: List[Dict]) -> Dict:
        """Categorize validated variants by risk level"""
        
        categories = {
            'low': [],
            'medium': [],
            'high': []
        }
        
        for v in variants:
            risk = v.get('risk_level', 'HIGH').lower()
            if risk in categories:
                categories[risk].append(v)
        
        return categories
    
    def _select_final_variants(self, variants: List[Dict], top_k: int = 15) -> List[Dict]:
        """Select final diverse variants"""
        
        if not variants:
            return []
        
        # Sort by consensus score
        variants.sort(key=lambda x: x.get('consensus_score', 0), reverse=True)
        
        # Ensure diversity across tracks
        final = []
        tracks_seen = {}
        
        for v in variants:
            track = v.get('track', 'unknown')
            
            # Limit per track to ensure diversity
            if track not in tracks_seen:
                tracks_seen[track] = 0
            
            if tracks_seen[track] < top_k // 3 + 1:
                final.append(v)
                tracks_seen[track] += 1
            
            if len(final) >= top_k:
                break
        
        return final
    
    def _generate_outputs(self, variants: List[Dict], metadata: Dict) -> Dict:
        """Generate all output files"""
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        outputs = {}
        
        # 1. FASTA file
        fasta_path = os.path.join(self.results_dir, f"d1_variants_{timestamp}.fasta")
        with open(fasta_path, 'w') as f:
            for i, v in enumerate(variants):
                validation = v.get('validation', {})
                f.write(f">{v.get('id', f'VAR_{i:04d}')} ")
                f.write(f"score={v.get('consensus_score', 0):.2f} ")
                f.write(f"delta_tm={validation.get('delta_tm', 0):+.1f} ")
                f.write(f"risk={v.get('risk_level', 'UNKNOWN')}\n")
                f.write(f"{v['sequence']}\n")
        
        outputs['fasta'] = fasta_path
        print(f"  ✓ FASTA: {fasta_path}")
        
        # 2. Detailed CSV
        csv_path = os.path.join(self.results_dir, f"d1_variants_{timestamp}.csv")
        df_data = []
        
        for v in variants:
            validation = v.get('validation', {})
            df_data.append({
                'ID': v.get('id'),
                'Track': v.get('track'),
                'Method': v.get('method'),
                'Sequence_Length': len(v['sequence']),
                'Num_Mutations': len(v.get('mutations', [])),
                'Consensus_Score': v.get('consensus_score', 0),
                'Predicted_Tm': v.get('predicted_tm', 0),
                'Delta_Tm': validation.get('delta_tm', 0),
                'pLDDT': validation.get('plddt', 0),
                'TM_Score': validation.get('tm_score', 0),
                'Aggregation_Score': validation.get('aggregation_score', 0),
                'Risk_Level': v.get('risk_level'),
                'Validation_Passes': validation.get('pass_rate', 0)
            })
        
        df = pd.DataFrame(df_data)
        df.to_csv(csv_path, index=False)
        outputs['csv'] = csv_path
        print(f"  ✓ CSV: {csv_path}")
        
        # 3. JSON metadata
        json_path = os.path.join(self.results_dir, f"d1_metadata_{timestamp}.json")
        with open(json_path, 'w') as f:
            json.dump(metadata, f, indent=2, default=str)
        
        outputs['metadata'] = json_path
        print(f"  ✓ Metadata: {json_path}")
        
        # 4. Validation reports for top 5
        for i, v in enumerate(variants[:5]):
            if 'validation' in v:
                report = self.validator.generate_report(v)
                report_path = os.path.join(
                    self.results_dir,
                    f"validation_report_{v.get('id', f'VAR_{i:04d}')}_{timestamp}.md"
                )
                with open(report_path, 'w') as f:
                    f.write(report)
                
                if i == 0:
                    print(f"  ✓ Generated validation reports for top 5 variants")
        
        return outputs
    
    def _print_summary(self, variants: List[Dict], results: Dict):
        """Print final summary"""
        
        print(f"""
╔════════════════════════════════════════════════════╗
║              PIPELINE COMPLETE                     ║
╚════════════════════════════════════════════════════╝

📊 RESULTS SUMMARY:
  • Total variants generated: {results['stages']['generation']['total_generated']}
  • Passed screening: {results['stages'].get('screening', {}).get('passed', 0)}
  • Validated: {results['stages'].get('validation', {}).get('total_validated', 0)}
  • Final selection: {len(variants)}

🏆 TOP VARIANTS:
""")
        
        for i, v in enumerate(variants[:5], 1):
            validation = v.get('validation', {})
            print(f"  {i}. {v.get('id', 'Unknown')}")
            print(f"     Track: {v.get('track', 'unknown')}")
            print(f"     Method: {v.get('method', 'unknown')}")
            print(f"     Consensus Score: {v.get('consensus_score', 0):.2f}")
            print(f"     ΔTm: {validation.get('delta_tm', 0):+.1f}°C")
            print(f"     Risk: {v.get('risk_level', 'UNKNOWN')}")
            print()
        
        print(f"""
💰 API USAGE:
  • Total API calls: {self.client.api_calls}
  • Estimated cost: ${self.client.api_costs:.2f}

📁 OUTPUT FILES:
  • Results saved to: {self.results_dir}
  • Files generated: {len(results.get('output_files', {}))}

⏱️ RUNTIME: {results.get('runtime_seconds', 0)/60:.1f} minutes
""")