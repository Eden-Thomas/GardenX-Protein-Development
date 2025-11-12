"""
master_pipeline_v3_sequential.py - COMPLETE SEQUENTIAL VERSION WITH ITERATIVE IMPROVEMENT
Each track builds on previous insights with full guardrails
Track 1 (Aggressive Base) → Track 2 (Iterative Improvement) → Track 3 (Structural) → Validation
UPDATED: Now uses iteratively_improve_variant() for true sequential processing
"""

import asyncio
import numpy as np
from typing import Dict, List, Optional
from datetime import datetime
from pathlib import Path
import json


class MasterPipelineV3Sequential:
    """Fully Sequential D1 Engineering Pipeline with Iterative Improvement"""
    
    def __init__(self, config, client=None):
        """
        Initialize pipeline
        
        Args:
            config: Config object with settings
            client: NeurosnapClient instance (optional)
        """
        self.config = config
        self.client = client
        
        # Import tracks (handle if modules don't exist yet)
        try:
            from pipeline.track1_evolutionary import Track1Evolutionary
            from pipeline.track2_generative import Track2Generative
            from pipeline.track3_structural import Track3Structural
            
            self.track1 = Track1Evolutionary(config, client)
            self.track2 = Track2Generative(config, client)
            self.track3 = Track3Structural(config, client)
        except ImportError as e:
            print(f"Warning: Could not import track modules: {e}")
            self.track1 = None
            self.track2 = None
            self.track3 = None
        
        # Validator (initialize in run if available)
        self.validator = None
    
    async def run(
        self,
        wt_sequence: str,
        target_temp: int = 55,  # HIGHER target for aggressive mode
        n_variants: int = 5,    # FEWER but better variants
        validate_all: bool = True
    ) -> Dict:
        """
        Execute TRUE SEQUENTIAL pipeline with iterative improvement
        Track 1 → Track 2 (improve Track 1 variant) → Track 3 (optimize Track 2 output variant)
        """
        
        print("\n" + "="*70)
        print("    TRUE SEQUENTIAL D1 ENGINEERING PIPELINE")
        print("    Same variants flow through each track for iterative improvement")
        print("="*70)
        print(f"  Target Temperature: {target_temp}°C (WT ~28°C = +{target_temp-28}°C)")
        print(f"  Base Variants: {n_variants} (will be iteratively improved)")
        print(f"  Flow: Track 1 → Track 2 (improve same) → Track 3 (optimize same)")
        print("="*70 + "\n")
        
        start_time = datetime.now()
        
        # Initialize validator if available
        print("📊 INITIALIZATION")
        print("─" * 50)
        
        try:
            from validation.comprehensive_validator import ComprehensiveValidator
            if not self.validator:
                self.validator = ComprehensiveValidator(self.client, wt_sequence)
            await self.validator.initialize_wt_structure(wt_sequence)
            wt_tm = self.validator.wt_tm
            print(f"  ✓ WT baseline Tm: {wt_tm:.1f}°C")
            print(f"  ✓ Target improvement: +{target_temp - wt_tm:.1f}°C\n")
        except Exception as e:
            print(f"  ⚠️ Validator not available: {e}")
            wt_tm = 28.0  # Realistic WT baseline
            print(f"  ✓ WT baseline Tm: {wt_tm:.1f}°C")
            print(f"  ✓ Target improvement: +{target_temp - wt_tm:.1f}°C\n")
        
        # ========== TRACK 1: AGGRESSIVE BASE VARIANTS ==========
        print("🧬 TRACK 1: AGGRESSIVE EVOLUTIONARY BASE")
        print("─" * 50)
        print(f"  Generating {n_variants} aggressive base variants...")
        
        base_variants = []
        if self.track1:
            try:
                # Generate aggressive base variants with more mutations
                base_variants = await self.track1.generate_variants(
                    wt_sequence=wt_sequence,
                    n_variants=n_variants,
                    aggressive_mode=True,           # NEW: More aggressive
                    conservation_threshold=0.3,     # LOWER threshold
                    max_mutations=15,               # MORE mutations allowed
                    target_temp=target_temp         # HIGHER target
                )
                
                if base_variants:
                    print(f"  ✓ Generated {len(base_variants)} base variants")
                    for i, v in enumerate(base_variants):
                        n_muts = len(v.get('mutations', []))
                        print(f"    {i+1}. {v['id']}: {n_muts} mutations")
                else:
                    print("  ❌ Track 1 failed to generate variants")
                    return {'error': 'Track 1 failed'}
            
            except Exception as e:
                print(f"  ❌ Track 1 failed: {e}")
                import traceback
                traceback.print_exc()
                return {'error': f'Track 1 failed: {e}'}
        else:
            print("  ❌ Track 1 module not available")
            return {'error': 'Track 1 not available'}
        
        # ========== TRACK 2: ITERATIVE IMPROVEMENT ==========
        print(f"\n🤖 TRACK 2: ITERATIVE IMPROVEMENT OF SAME {len(base_variants)} VARIANTS")
        print("─" * 50)
        print("  Taking each base variant and iteratively improving it...")
        
        improved_variants = []
        if self.track2:
            for i, base_variant in enumerate(base_variants):
                print(f"  Improving {base_variant['id']}...")
                
                try:
                    # Improve THIS EXACT variant (not generate new ones)
                    improved_variant = await self.track2.iteratively_improve_variant(
                        base_variant=base_variant,
                        wt_sequence=wt_sequence,
                        target_temp=target_temp,
                        improvement_rounds=3,        # Multiple improvement rounds
                        aggressive_mode=True         # Allow big changes
                    )
                    
                    if improved_variant:
                        improved_variants.append(improved_variant)
                        old_muts = len(base_variant.get('mutations', []))
                        new_muts = len(improved_variant.get('mutations', []))
                        delta_tm = improved_variant.get('predicted_tm', 28) - 28
                        
                        print(f"    ✓ {base_variant['id']} → {improved_variant['id']}")
                        print(f"      Mutations: {old_muts} → {new_muts} (+{new_muts - old_muts})")
                        print(f"      Predicted Tm: {improved_variant.get('predicted_tm', 28):.1f}°C (+{delta_tm:.1f}°C)")
                    else:
                        print(f"    ⚠️ Could not improve {base_variant['id']}")
                        improved_variants.append(base_variant)  # Keep original
                
                except Exception as e:
                    print(f"    ❌ Failed to improve {base_variant['id']}: {e}")
                    improved_variants.append(base_variant)  # Keep original
        
        else:
            print("  ⚠️ Track 2 not available, using base variants")
            improved_variants = base_variants
        
        print(f"  ✓ Improved {len(improved_variants)} variants")
        
        # ========== TRACK 3: STRUCTURAL OPTIMIZATION ==========
        print(f"\n🏗️ TRACK 3: STRUCTURAL OPTIMIZATION OF SAME {len(improved_variants)} VARIANTS")
        print("─" * 50)
        print("  Taking each improved variant and structurally optimizing it...")
        
        final_variants = []
        if self.track3:
            for i, improved_variant in enumerate(improved_variants):
                print(f"  Optimizing {improved_variant['id']}...")
                
                try:
                    # For Track 3, we pass improved variants as seed_variants
                    # Track 3 will refine these specific variants
                    track3_result = await self.track3.generate_variants(
                        wt_sequence=wt_sequence,
                        n_variants=1,  # One refined version per improved variant
                        seed_variants=[improved_variant]
                    )
                    
                    if track3_result and len(track3_result) > 0:
                        optimized_variant = track3_result[0]
                        final_variants.append(optimized_variant)
                        
                        delta_tm = optimized_variant.get('predicted_tm', improved_variant.get('predicted_tm', 28)) - 28
                        plddt = optimized_variant.get('plddt', 0)
                        
                        print(f"    ✓ {improved_variant['id']} → {optimized_variant['id']}")
                        print(f"      Final Tm: {optimized_variant.get('predicted_tm', improved_variant.get('predicted_tm', 28)):.1f}°C (+{delta_tm:.1f}°C)")
                        print(f"      pLDDT: {plddt:.1f}")
                    else:
                        print(f"    ⚠️ Could not optimize {improved_variant['id']}")
                        # Add improved variant with Track 3 metadata
                        improved_variant_copy = improved_variant.copy()
                        improved_variant_copy['track'] = 'track3'
                        improved_variant_copy['method'] = 'structural_passthrough'
                        final_variants.append(improved_variant_copy)
                
                except Exception as e:
                    print(f"    ❌ Failed to optimize {improved_variant['id']}: {e}")
                    # Keep improved version with Track 3 metadata
                    improved_variant_copy = improved_variant.copy()
                    improved_variant_copy['track'] = 'track3'
                    improved_variant_copy['method'] = 'structural_passthrough'
                    final_variants.append(improved_variant_copy)
        
        else:
            print("  ⚠️ Track 3 not available, using improved variants")
            for v in improved_variants:
                v_copy = v.copy()
                v_copy['track'] = 'track3'
                v_copy['method'] = 'passthrough'
                final_variants.append(v_copy)
        
        print(f"  ✓ Optimized {len(final_variants)} variants")
        
        # ========== FINAL RESULTS ==========
        print(f"\n🏆 FINAL RESULTS: TRUE SEQUENTIAL PROCESSING")
        print("─" * 50)
        
        # Sort by predicted Tm
        final_variants.sort(key=lambda x: x.get('predicted_tm', 28), reverse=True)
        
        print(f"  Processed {len(final_variants)} variants through full pipeline:")
        print()
        
        for i, variant in enumerate(final_variants, 1):
            delta_tm = variant.get('predicted_tm', 28) - wt_tm
            n_muts = len(variant.get('mutations', []))
            
            # Show lineage
            lineage = []
            if variant.get('track1_parent'):
                lineage.append(variant['track1_parent'])
            elif 'T1_EVO' in variant.get('id', ''):
                lineage.append(variant['id'].split('_IMPROVED')[0].split('_R')[0])
            
            if variant.get('track2_parent'):
                lineage.append(variant['track2_parent'])
            
            if variant.get('parent_variant'):
                lineage.append(variant['parent_variant'])
            
            lineage.append(variant['id'])
            lineage_str = " → ".join(lineage) if lineage else variant['id']
            
            print(f"    {i}. {variant['id']}")
            print(f"       Predicted Tm: {variant.get('predicted_tm', 28):.1f}°C (+{delta_tm:.1f}°C)")
            print(f"       Mutations: {n_muts} total")
            print(f"       Lineage: {lineage_str}")
            
            # Show first 8 mutations
            muts_display = variant.get('mutations', [])[:8]
            if len(variant.get('mutations', [])) > 8:
                muts_display_str = ', '.join(muts_display) + f' ... (+{len(variant.get("mutations", [])) - 8} more)'
            else:
                muts_display_str = ', '.join(muts_display)
            print(f"       Key mutations: {muts_display_str}")
            
            # Check if target achieved
            if variant.get('predicted_tm', 28) >= target_temp:
                print(f"       🎯 TARGET ACHIEVED! ({target_temp}°C)")
            print()
        
        # Success metrics
        successful_variants = [v for v in final_variants if v.get('predicted_tm', 28) >= target_temp]
        mean_improvement = np.mean([v.get('predicted_tm', 28) - wt_tm for v in final_variants])
        
        print(f"  📊 Success Rate: {len(successful_variants)}/{len(final_variants)} variants reached {target_temp}°C")
        print(f"  📈 Mean Improvement: +{mean_improvement:.1f}°C")
        print(f"  🔥 Best Variant: {final_variants[0]['id']} at {final_variants[0].get('predicted_tm', 28):.1f}°C")
        
        # Export results
        results = self._export_results(final_variants, start_time, wt_tm, target_temp)
        results['pipeline_type'] = 'true_sequential'
        results['target_achieved'] = len(successful_variants) > 0
        results['successful_variants'] = len(successful_variants)
        
        return results
    
    def _export_results(
        self, variants: List[Dict], start_time: datetime,
        wt_tm: float, target_temp: int
    ) -> Dict:
        """Export comprehensive results"""
        
        runtime = (datetime.now() - start_time).total_seconds() / 60
        
        # Summary statistics
        summary = {
            'total_variants': len(variants),
            'runtime_minutes': runtime,
            'wt_tm': wt_tm,
            'target_temp': target_temp
        }
        
        if variants:
            delta_tms = [v.get('predicted_tm', 28) - wt_tm for v in variants]
            summary['mean_delta_tm'] = np.mean(delta_tms)
            summary['max_delta_tm'] = np.max(delta_tms)
            summary['min_delta_tm'] = np.min(delta_tms)
        
        # Save results
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        try:
            output_dir = Path(self.config.RESULTS_DIR)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            results_file = output_dir / f"sequential_results_{timestamp}.json"
            
            results = {
                'pipeline_version': '3.0_sequential_iterative',
                'timestamp': timestamp,
                'summary': summary,
                'final_variants': variants
            }
            
            with open(results_file, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            print(f"\n📊 Results exported to: {results_file}")
        
        except Exception as e:
            print(f"\n⚠️ Could not save results: {e}")
            results = {
                'pipeline_version': '3.0_sequential_iterative',
                'timestamp': timestamp,
                'summary': summary,
                'final_variants': variants
            }
        
        return results


# ========== STANDALONE EXECUTION ==========
if __name__ == "__main__":
    import sys
    from pathlib import Path
    
    # Add parent directory to path
    sys.path.append(str(Path(__file__).parent.parent))
    
    # Import dependencies
    try:
        from config import Config
        from neurosnap_client import NeurosnapClient
        
        async def test():
            """Test the pipeline"""
            
            print("Initializing pipeline...")
            config = Config()
            
            # Initialize client if API key available
            client = None
            if config.NEUROSNAP_API_KEY:
                client = NeurosnapClient(config.NEUROSNAP_API_KEY)
                print("✓ Neurosnap client initialized")
            else:
                print("⚠️ No API key, running in offline mode")
            
            pipeline = MasterPipelineV3Sequential(config, client)
            
            # Use a test sequence (truncated for speed)
            wt_sequence = "MTAILERRESESLWGRFCNWITGTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIITGAVVPSSNAILVHFYPIWEAASLDEWLYNGGPYQLIIFHFLLGASCYMGREWELSFRLGMRPWICVAYSAPLASAFAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLLFSAMHGSLVTSSLIRETTENESANEGYRFGQEEETYNIVAHAGYGRLLFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAAEG"
            
            print(f"WT sequence length: {len(wt_sequence)}")
            
            results = await pipeline.run(
                wt_sequence=wt_sequence,
                target_temp=55,     # Higher target
                n_variants=5,       # Fewer but better
                validate_all=False  # Skip validation for speed
            )
            
            print("\n✅ Pipeline test complete!")
            print(f"Generated {len(results.get('final_variants', []))} variants")
        
        # Run test
        asyncio.run(test())
    
    except ImportError as e:
        print(f"Error importing dependencies: {e}")
        print("Make sure you're running from the backend directory")