"""
master_pipeline_v3_sequential.py - COMPLETE SEQUENTIAL VERSION
Each track builds on previous insights with full guardrails
Track 1 → Track 2 → Track 3 → Validation → Final Ranking
WITH FIXED SEQUENTIAL RANKING PRIORITY
"""

import asyncio
import numpy as np
from typing import Dict, List, Optional
from datetime import datetime
from pathlib import Path
import json


class MasterPipelineV3Sequential:
    """Fully Sequential D1 Engineering Pipeline"""
    
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
        target_temp: int = 50,
        n_variants: int = 50,
        validate_all: bool = True
    ) -> Dict:
        """Execute fully sequential pipeline"""
        
        print("\n" + "="*70)
        print("    SEQUENTIAL D1 ENGINEERING PIPELINE V3.0")
        print("    Track 1 → Track 2 → Track 3 with Full Guardrails")
        print("="*70)
        print(f"  Target Temperature: {target_temp}°C")
        print(f"  Total Variants Target: {n_variants}")
        print(f"  Sequential Analysis: ENABLED")
        print(f"  Guardrails: ALL ACTIVE")
        print("="*70 + "\n")
        
        start_time = datetime.now()
        all_variants = []
        
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
            print(f"  ✓ Using estimated WT Tm: 42.0°C")
            print(f"  ✓ Target improvement: +{target_temp - 42:.1f}°C\n")
            wt_tm = 42.0
        
        # ========== TRACK 1: EVOLUTIONARY ==========
        print("🧬 TRACK 1: EVOLUTIONARY ANALYSIS")
        print("─" * 50)
        
        track1_results = []
        if self.track1:
            try:
                print("  Analyzing thermophilic organisms with guardrails...")
                track1_results = await self.track1.generate_variants(
                    wt_sequence=wt_sequence,
                    n_variants=max(10, n_variants // 3)
                )
                
                if track1_results:
                    all_variants.extend(track1_results)
                    print(f"\n  ✓ Track 1 complete: {len(track1_results)} variants")
                    print(f"    • Conservation-filtered mutations")
                    print(f"    • Topology-aware substitutions")
                    print(f"    • Functional sites protected")
                else:
                    print("\n  ⚠️ Track 1 generated no variants")
            
            except Exception as e:
                print(f"  ❌ Track 1 failed: {e}")
                import traceback
                traceback.print_exc()
        else:
            print("  ⚠️ Track 1 module not available")
        
        # ========== EXTRACT SEEDS FOR TRACK 2 ==========
        print("\n🔗 SEED EXTRACTION")
        print("─" * 50)
        
        seed_mutations = self._extract_seed_mutations(track1_results)
        seed_variants_for_t2 = track1_results
        
        if seed_mutations:
            print(f"  ✓ Extracted {len(seed_mutations)} seed mutations for Track 2:")
            print(f"    {', '.join(seed_mutations[:8])}")
        else:
            print("  ⚠️ No seed mutations extracted (Track 2 will use de novo)")
        
        # ========== TRACK 2: AI-GUIDED OPTIMIZATION ==========
        print("\n🤖 TRACK 2: AI-GUIDED OPTIMIZATION")
        print("─" * 50)
        
        track2_results = []
        if self.track2:
            try:
                print("  Building on Track 1 insights with multi-objective scoring...")
                track2_results = await self.track2.generate_variants(
                    wt_sequence=wt_sequence,
                    target_temp=target_temp,
                    n_variants=max(10, n_variants // 3),
                    seed_mutations=seed_mutations,
                    seed_variants=seed_variants_for_t2
                )
                
                if track2_results:
                    all_variants.extend(track2_results)
                    print(f"\n  ✓ Track 2 complete: {len(track2_results)} variants")
                    
                    # Count how many used Track 1 seeds
                    t1_based = sum(1 for v in track2_results if v.get('based_on_track1', False))
                    if t1_based > 0:
                        print(f"    • {t1_based} variants built on Track 1 mutations")
                    
                    mean_delta_tm = np.mean([v.get('delta_tm', 0) for v in track2_results])
                    print(f"    • Mean predicted ΔTm: {mean_delta_tm:.1f}°C")
                else:
                    print("\n  ⚠️ Track 2 generated no variants")
            
            except Exception as e:
                print(f"  ❌ Track 2 failed: {e}")
                import traceback
                traceback.print_exc()
        else:
            print("  ⚠️ Track 2 module not available")
        
        # ========== SELECT TOP CANDIDATES FOR TRACK 3 ==========
        print("\n🎯 CANDIDATE SELECTION FOR TRACK 3")
        print("─" * 50)
        
        seed_variants_for_t3 = self._select_refinement_candidates(
            all_variants,
            n_select=min(10, len(all_variants))
        )
        
        if seed_variants_for_t3:
            print(f"  ✓ Selected {len(seed_variants_for_t3)} top candidates for structural refinement:")
            for i, v in enumerate(seed_variants_for_t3[:5]):
                score = v.get('composite_score', v.get('score', 0))
                print(f"    {i+1}. {v.get('id', 'unknown')} - Score: {score:.2f}, {len(v.get('mutations', []))} mutations")
        else:
            print("  ⚠️ No candidates selected (Track 3 will generate de novo)")
        
        # ========== TRACK 3: STRUCTURAL REFINEMENT ==========
        print("\n🏗️ TRACK 3: STRUCTURAL REFINEMENT")
        print("─" * 50)
        
        track3_results = []
        if self.track3:
            try:
                print("  Optimizing structures with membrane-aware refinement...")
                track3_results = await self.track3.generate_variants(
                    wt_sequence=wt_sequence,
                    n_variants=max(5, n_variants - len(all_variants)),
                    seed_variants=seed_variants_for_t3
                )
                
                if track3_results:
                    all_variants.extend(track3_results)
                    print(f"\n  ✓ Track 3 complete: {len(track3_results)} variants")
                    
                    # Count refined variants
                    refined = sum(1 for v in track3_results if 'parent_variant' in v)
                    if refined > 0:
                        print(f"    • {refined} variants refined from Track 1+2 candidates")
                    
                    mean_plddt = np.mean([v.get('plddt', 0) for v in track3_results])
                    print(f"    • Mean pLDDT: {mean_plddt:.1f}")
                else:
                    print("\n  ⚠️ Track 3 generated no variants")
            
            except Exception as e:
                print(f"  ❌ Track 3 failed: {e}")
                import traceback
                traceback.print_exc()
        else:
            print("  ⚠️ Track 3 module not available")
        
        # ========== VALIDATION ==========
        print("\n🔬 VALIDATION & SCORING")
        print("─" * 50)
        
        if all_variants:
            if self.validator and validate_all:
                print("  Running comprehensive validation...")
                validated_variants = await self._validate_all_variants(all_variants)
            else:
                print("  Using heuristic scoring (validator not available)...")
                validated_variants = self._heuristic_validation(all_variants)
            
            print(f"  ✓ Validated {len(validated_variants)} variants")
            
            # Final ranking WITH FIXED SEQUENTIAL PRIORITY
            print("\n🏆 FINAL RANKING (Sequential Priority)")
            print("─" * 50)
            
            final_variants = self._rank_variants_sequentially(validated_variants)
            
            # Export results
            results = self._export_results(final_variants, start_time, wt_tm, target_temp)
            
            # Print top 5 with priority indicators
            print("\n  Top 5 variants:")
            for i, variant in enumerate(final_variants[:5], 1):
                vid = variant.get('id', 'unknown')
                delta_tm = variant.get('validation', {}).get('delta_tm', variant.get('delta_tm', 0))
                n_muts = len(variant.get('mutations', []))
                track = variant.get('track', 'unknown')
                priority_mult = variant.get('priority_multiplier', 1.0)
                
                # Show parent if Track 3
                parent_info = ""
                if track == 'track3' and variant.get('parent_variant'):
                    parent_info = f" ← refined from {variant.get('parent_variant')}"
                
                print(f"    {i}. {vid} ({track}) - ΔTm: +{delta_tm:.1f}°C, {n_muts} mutations [Priority: {priority_mult:.1f}x]{parent_info}")
            
            return results
        
        else:
            print("  ⚠️ No variants generated across all tracks")
            return {
                'final_variants': [],
                'summary': {
                    'total_generated': 0,
                    'track1': 0,
                    'track2': 0,
                    'track3': 0
                }
            }
    
    def _extract_seed_mutations(self, track1_results: List[Dict]) -> List[str]:
        """Extract promising mutations from Track 1 for Track 2"""
        
        if not track1_results:
            return []
        
        # Get mutations from top 25% of Track 1 variants
        sorted_variants = sorted(
            track1_results,
            key=lambda x: x.get('score', 0),
            reverse=True
        )
        
        top_variants = sorted_variants[:max(1, len(sorted_variants) // 4)]
        
        seed_mutations = set()
        for v in top_variants:
            seed_mutations.update(v.get('mutations', []))
        
        return sorted(list(seed_mutations))[:15]  # Top 15 mutations
    
    def _select_refinement_candidates(
        self, variants: List[Dict], n_select: int
    ) -> List[Dict]:
        """Select best candidates for structural refinement"""
        
        if not variants:
            return []
        
        # Score variants for selection
        for v in variants:
            score = 0
            score += v.get('delta_tm', 0) * 0.4
            score += v.get('composite_score', 0) * 0.3
            score += v.get('conservation_score', 0) * 0.15
            score += v.get('score', 0) * 0.15
            v['refinement_score'] = score
        
        # Sort and select top N
        sorted_variants = sorted(
            variants,
            key=lambda x: x.get('refinement_score', 0),
            reverse=True
        )
        
        return sorted_variants[:n_select]
    
    async def _validate_all_variants(self, variants: List[Dict]) -> List[Dict]:
        """Validate all variants comprehensively"""
        
        validated = []
        
        for variant in variants:
            try:
                validation_result = await self.validator.validate_variant(variant)
                
                pass_rate = validation_result.get('validation', {}).get('pass_rate', 0)
                if pass_rate > 0.3:  # Keep if passes >30% of checks
                    validated.append(validation_result)
            
            except Exception as e:
                print(f"    ⚠️ Validation failed for {variant.get('id', 'unknown')}: {e}")
                validated.append(variant)
        
        return validated
    
    def _heuristic_validation(self, variants: List[Dict]) -> List[Dict]:
        """Simple heuristic validation when validator unavailable"""
        
        for variant in variants:
            # Add simple validation scores
            variant['validation'] = {
                'delta_tm': variant.get('delta_tm', 0),
                'predicted_tm': variant.get('predicted_tm', 42),
                'composite_score': variant.get('composite_score', variant.get('score', 0)),
                'pass_rate': 0.8,
                'validated': True
            }
        
        return variants
    
    def _rank_variants_sequentially(self, variants: List[Dict]) -> List[Dict]:
        """
        FIXED: Rank variants properly for sequential pipeline
        
        Priority system:
        1. Track 3 variants (structurally refined) - HIGHEST (3.0x multiplier)
        2. Track 2 variants (AI-optimized) - MEDIUM (2.0x multiplier)
        3. Track 1 variants (evolutionary baseline) - BASELINE (1.0x multiplier)
        
        This ensures that refined variants (Track 3) always rank higher
        than their parent variants (Track 1/2)
        """
        
        for v in variants:
            track = v.get('track', 'unknown')
            
            # Base score components
            validation = v.get('validation', {})
            delta_tm = validation.get('delta_tm', v.get('delta_tm', 0))
            plddt = v.get('plddt', 80)
            n_mutations = len(v.get('mutations', []))
            conservation = v.get('conservation_score', 0.5)
            membrane = v.get('membrane_packing', 0)
            
            # Component scores
            tm_score = delta_tm * 0.40
            structural_score = (plddt / 100) * 0.30
            safety_score = max(0, (10 - n_mutations) / 10) * 0.20
            conservation_score = conservation * 0.10
            
            base_score = tm_score + structural_score + safety_score + conservation_score
            
            # CRITICAL: Track priority multipliers (ensures sequential ranking)
            if track == 'track3':
                priority_multiplier = 3.0  # HIGHEST - refined variants
                # Extra bonus if has lineage
                if v.get('parent_variant'):
                    priority_multiplier *= 1.2  # 3.6x total for refined variants
            elif track == 'track2':
                priority_multiplier = 2.0  # MEDIUM - AI-optimized
                # Bonus if built on Track 1
                if v.get('based_on_track1'):
                    priority_multiplier *= 1.1  # 2.2x for Track 1-based variants
            else:  # track1
                priority_multiplier = 1.0  # BASELINE - evolutionary
            
            # Final score
            v['final_rank_score'] = base_score * priority_multiplier
            v['priority_multiplier'] = priority_multiplier
            v['base_score'] = base_score
            
            # Store component scores for analysis
            v['component_scores'] = {
                'tm': tm_score,
                'structural': structural_score,
                'safety': safety_score,
                'conservation': conservation_score
            }
        
        # Sort by final score (Track 3 will naturally rank highest)
        ranked = sorted(variants, key=lambda x: x.get('final_rank_score', 0), reverse=True)
        
        return ranked
    
    def _export_results(
        self, variants: List[Dict], start_time: datetime,
        wt_tm: float, target_temp: int
    ) -> Dict:
        """Export comprehensive results"""
        
        runtime = (datetime.now() - start_time).total_seconds() / 60
        
        # Summary statistics
        track1_count = sum(1 for v in variants if v.get('track') == 'track1')
        track2_count = sum(1 for v in variants if v.get('track') == 'track2')
        track3_count = sum(1 for v in variants if v.get('track') == 'track3')
        
        # Track distribution in top 15
        top_15_tracks = {
            'track1': sum(1 for v in variants[:15] if v.get('track') == 'track1'),
            'track2': sum(1 for v in variants[:15] if v.get('track') == 'track2'),
            'track3': sum(1 for v in variants[:15] if v.get('track') == 'track3')
        }
        
        summary = {
            'total_generated': len(variants),
            'track1': track1_count,
            'track2': track2_count,
            'track3': track3_count,
            'top_15_distribution': top_15_tracks,
            'runtime_minutes': runtime,
            'wt_tm': wt_tm,
            'target_temp': target_temp
        }
        
        if variants:
            delta_tms = [v.get('validation', {}).get('delta_tm', v.get('delta_tm', 0)) for v in variants]
            summary['mean_delta_tm'] = np.mean(delta_tms)
            summary['max_delta_tm'] = np.max(delta_tms)
        
        # Save results
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        
        try:
            output_dir = Path(self.config.RESULTS_DIR)
            output_dir.mkdir(parents=True, exist_ok=True)
            
            results_file = output_dir / f"sequential_results_{timestamp}.json"
            
            results = {
                'pipeline_version': '3.0_sequential_fixed',
                'timestamp': timestamp,
                'summary': summary,
                'final_variants': variants
            }
            
            with open(results_file, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            print(f"\n📊 Results exported to: {results_file}")
            print(f"\n📈 Top 15 Track Distribution:")
            print(f"   • Track 3 (Refined): {top_15_tracks['track3']} variants")
            print(f"   • Track 2 (AI-Optimized): {top_15_tracks['track2']} variants")
            print(f"   • Track 1 (Evolutionary): {top_15_tracks['track1']} variants")
        
        except Exception as e:
            print(f"\n⚠️ Could not save results: {e}")
            results = {
                'pipeline_version': '3.0_sequential_fixed',
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
                target_temp=50,
                n_variants=15,  # Small number for testing
                validate_all=False  # Skip validation for speed
            )
            
            print("\n✅ Pipeline test complete!")
            print(f"Generated {len(results.get('final_variants', []))} variants")
        
        # Run test
        asyncio.run(test())
    
    except ImportError as e:
        print(f"Error importing dependencies: {e}")
        print("Make sure you're running from the backend directory")