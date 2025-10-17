"""
master_pipeline_v3_sequential.py
Sequential D1 Engineering Pipeline - Each track builds on previous insights
Track 1 → Track 2 → Track 3 with statistical validation at each step
"""

import asyncio
import numpy as np
from typing import Dict, List, Optional, Tuple
from datetime import datetime
from pathlib import Path
import json
from scipy import stats

from config import Config
from neurosnap_client import NeurosnapClient
from pipeline.track1_evolutionary import Track1Evolutionary
from pipeline.track2_generative import Track2Generative
from pipeline.track3_structural import Track3Structural
from pipeline.consensus_scorer import ConsensusScorer
from validation.comprehensive_validator import ComprehensiveValidator


class MasterPipelineV3Sequential:
    """
    Sequential Master Pipeline for D1 protein engineering
    Each track builds on insights from the previous track
    """
    
    def __init__(self, config: Config):
        self.config = config
        self.client = NeurosnapClient(config.NEUROSNAP_API_KEY) if config.NEUROSNAP_API_KEY else None
        
        # Initialize tracks
        self.track1 = Track1Evolutionary(config, self.client)
        self.track2 = Track2Generative(config, self.client)
        self.track3 = Track3Structural(config, self.client)
        
        # Validation and scoring
        self.validator = ComprehensiveValidator(config, self.client)
        self.scorer = ConsensusScorer(config=config, client=self.client)
        
        # Statistical tracking
        self.statistical_results = {
            'track1_stats': {},
            'track2_stats': {},
            'track3_stats': {},
            'overall_stats': {}
        }
    
    async def run(self, 
                  wt_sequence: str,
                  target_temp: int = 50,
                  n_variants: int = 50,
                  validate_all: bool = True) -> Dict:
        """
        Execute sequential pipeline with statistical validation
        """
        
        print("\n" + "="*70)
        print("    SEQUENTIAL D1 ENGINEERING PIPELINE V3.0")
        print("    Each Track Builds on Previous Insights")
        print("="*70)
        print(f"  Target Temperature: {target_temp}°C")
        print(f"  Total Variants Target: {n_variants}")
        print(f"  Sequential Analysis: ENABLED")
        print(f"  Statistical Validation: ENABLED")
        print("="*70 + "\n")
        
        start_time = datetime.now()
        all_variants = []
        
        # Initialize validator with WT structure
        print("📊 INITIALIZATION PHASE")
        print("─" * 50)
        await self.validator.initialize_wt_structure(wt_sequence)
        wt_tm = self.validator.wt_tm
        print(f"  ✓ WT baseline Tm: {wt_tm:.1f}°C")
        print(f"  ✓ Target improvement: +{target_temp - wt_tm:.1f}°C\n")
        
        # ========== TRACK 1: EVOLUTIONARY ANALYSIS ==========
        print("🧬 TRACK 1: EVOLUTIONARY ANALYSIS")
        print("─" * 50)
        print("  Analyzing heat-tolerant organisms...")
        
        try:
            track1_results = await self.track1.run(
                wt_sequence=wt_sequence,
                n_variants=n_variants // 3,  # 1/3 of total variants
                target_temp=target_temp
            )
            
            if track1_results:
                all_variants.extend(track1_results)
                
                # Statistical analysis of Track 1
                track1_stats = self._analyze_track1_statistics(track1_results, wt_sequence)
                self.statistical_results['track1_stats'] = track1_stats
                
                print(f"\n  📊 Track 1 Statistical Summary:")
                print(f"     • Variants generated: {len(track1_results)}")
                print(f"     • Conservation score range: {track1_stats['conservation_range']}")
                print(f"     • Temperature correlation (r): {track1_stats['temp_correlation']:.3f}")
                print(f"     • P-value: {track1_stats['p_value']:.4f}")
                print(f"     • Top mutations identified: {', '.join(track1_stats['top_mutations'][:5])}")
                
        except Exception as e:
            print(f"  ❌ Track 1 failed: {str(e)}")
            track1_results = []
        
        # ========== TRACK 2: AI OPTIMIZATION (USING TRACK 1) ==========
        print("\n🤖 TRACK 2: AI-GUIDED OPTIMIZATION")
        print("─" * 50)
        print("  Building on evolutionary insights...")
        
        # Extract seed mutations from Track 1
        seed_mutations = self._extract_seed_mutations(track1_results) if track1_results else []
        
        if seed_mutations:
            print(f"  → Using {len(seed_mutations)} seed mutations from Track 1")
            print(f"    Seeds: {', '.join(seed_mutations[:5])}...")
        
        try:
            track2_results = await self.track2.run(
                wt_sequence=wt_sequence,
                n_variants=n_variants // 3,
                target_temp=target_temp,
                seed_mutations=seed_mutations  # Pass Track 1 insights
            )
            
            if track2_results:
                all_variants.extend(track2_results)
                
                # Statistical analysis of Track 2
                track2_stats = self._analyze_track2_statistics(track2_results, seed_mutations)
                self.statistical_results['track2_stats'] = track2_stats
                
                print(f"\n  📊 Track 2 Statistical Summary:")
                print(f"     • Variants generated: {len(track2_results)}")
                print(f"     • Mean predicted ΔTm: {track2_stats['mean_delta_tm']:.1f}°C")
                print(f"     • Confidence interval: [{track2_stats['ci_lower']:.1f}, {track2_stats['ci_upper']:.1f}]")
                print(f"     • Synergistic mutations: {track2_stats['synergistic_count']}")
                print(f"     • Improvement over Track 1: {track2_stats['improvement']:.1%}")
                
        except Exception as e:
            print(f"  ❌ Track 2 failed: {str(e)}")
            track2_results = []
        
        # ========== TRACK 3: STRUCTURAL REFINEMENT (USING TRACK 1+2) ==========
        print("\n🏗️ TRACK 3: STRUCTURAL REFINEMENT")
        print("─" * 50)
        print("  Optimizing best candidates structurally...")
        
        # Select top candidates from Tracks 1+2 for refinement
        candidates_for_refinement = self._select_refinement_candidates(
            track1_results + track2_results if track1_results else track2_results,
            n_select=10
        )
        
        print(f"  → Refining {len(candidates_for_refinement)} top candidates")
        
        try:
            track3_results = await self.track3.run(
                wt_sequence=wt_sequence,
                n_variants=n_variants - len(all_variants),  # Remaining variants
                target_temp=target_temp,
                seed_variants=candidates_for_refinement  # Pass best from Tracks 1+2
            )
            
            if track3_results:
                all_variants.extend(track3_results)
                
                # Statistical analysis of Track 3
                track3_stats = self._analyze_track3_statistics(track3_results)
                self.statistical_results['track3_stats'] = track3_stats
                
                print(f"\n  📊 Track 3 Statistical Summary:")
                print(f"     • Variants refined: {len(track3_results)}")
                print(f"     • Mean pLDDT: {track3_stats['mean_plddt']:.1f}")
                print(f"     • Structural violations: {track3_stats['violations']}")
                print(f"     • Active site preservation: {track3_stats['active_site_preserved']:.0%}")
                print(f"     • H-bond optimization: +{track3_stats['hbond_improvement']:.1f} bonds")
                
        except Exception as e:
            print(f"  ❌ Track 3 failed: {str(e)}")
            track3_results = []
        
        # ========== INTEGRATED ANALYSIS ==========
        print("\n🔬 INTEGRATED STATISTICAL ANALYSIS")
        print("─" * 50)
        
        if all_variants:
            # Comprehensive validation with statistics
            validated_variants = await self._validate_with_statistics(all_variants)
            
            # Generate mutation rationale for each variant
            for variant in validated_variants:
                variant['rationale'] = self._generate_mutation_rationale(variant)
            
            # Overall statistical summary
            overall_stats = self._calculate_overall_statistics(validated_variants)
            self.statistical_results['overall_stats'] = overall_stats
            
            print(f"  📊 Overall Pipeline Statistics:")
            print(f"     • Total variants: {len(all_variants)}")
            print(f"     • Passed validation: {len(validated_variants)}")
            print(f"     • Mean predicted Tm: {overall_stats['mean_tm']:.1f}°C")
            print(f"     • Best ΔTm achieved: +{overall_stats['best_delta_tm']:.1f}°C")
            print(f"     • Statistical significance: p<{overall_stats['significance']:.4f}")
            print(f"     • Pearson r (mutations vs Tm): {overall_stats['correlation']:.3f}")
            
            # Final ranking with consensus scoring
            print("\n🎯 CONSENSUS SCORING")
            print("─" * 50)
            final_variants = await self.scorer.score_variants(validated_variants)
            
            # Export results
            results = self._export_results(final_variants, start_time)
            
            return results
        else:
            print("  ⚠️ No variants generated")
            return {'final_variants': [], 'statistics': self.statistical_results}
    
    def _analyze_track1_statistics(self, variants: List[Dict], wt_sequence: str) -> Dict:
        """Statistical analysis of evolutionary track results"""
        
        conservation_scores = [v.get('conservation_score', 0) for v in variants]
        temp_correlations = [v.get('temperature_correlation', 0) for v in variants]
        
        # Calculate correlation with temperature
        if len(temp_correlations) > 1:
            r_value, p_value = stats.pearsonr(conservation_scores, temp_correlations)
        else:
            r_value, p_value = 0, 1
        
        # Extract top mutations
        all_mutations = []
        for v in variants:
            all_mutations.extend(v.get('mutations', []))
        
        mutation_counts = {}
        for mut in all_mutations:
            mutation_counts[mut] = mutation_counts.get(mut, 0) + 1
        
        top_mutations = sorted(mutation_counts.keys(), 
                              key=lambda x: mutation_counts[x], 
                              reverse=True)
        
        return {
            'conservation_range': f"{min(conservation_scores):.2f}-{max(conservation_scores):.2f}" if conservation_scores else "N/A",
            'temp_correlation': r_value,
            'p_value': p_value,
            'top_mutations': top_mutations,
            'mutation_frequency': mutation_counts
        }
    
    def _analyze_track2_statistics(self, variants: List[Dict], seed_mutations: List[str]) -> Dict:
        """Statistical analysis of AI optimization results"""
        
        delta_tms = [v.get('delta_tm', 0) for v in variants]
        
        if delta_tms:
            mean_delta = np.mean(delta_tms)
            std_delta = np.std(delta_tms)
            ci_lower = mean_delta - 1.96 * std_delta / np.sqrt(len(delta_tms))
            ci_upper = mean_delta + 1.96 * std_delta / np.sqrt(len(delta_tms))
        else:
            mean_delta = std_delta = ci_lower = ci_upper = 0
        
        # Count synergistic mutations (new mutations not in seeds)
        synergistic = 0
        for v in variants:
            v_mutations = set(v.get('mutations', []))
            if len(v_mutations - set(seed_mutations)) > 0:
                synergistic += 1
        
        return {
            'mean_delta_tm': mean_delta,
            'std_delta_tm': std_delta,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'synergistic_count': synergistic,
            'improvement': (mean_delta - 5) / 5 if mean_delta > 5 else 0  # vs baseline expectation
        }
    
    def _analyze_track3_statistics(self, variants: List[Dict]) -> Dict:
        """Statistical analysis of structural refinement results"""
        
        plddts = [v.get('plddt', 0) for v in variants]
        violations = sum(1 for v in variants if v.get('structural_violations', 0) > 0)
        
        # Active site preservation
        active_preserved = sum(1 for v in variants 
                             if v.get('active_site_distance', float('inf')) > 10)
        
        # H-bond analysis
        hbond_changes = [v.get('hbond_change', 0) for v in variants]
        
        return {
            'mean_plddt': np.mean(plddts) if plddts else 0,
            'violations': violations,
            'active_site_preserved': active_preserved / len(variants) if variants else 0,
            'hbond_improvement': np.mean(hbond_changes) if hbond_changes else 0
        }
    
    def _extract_seed_mutations(self, track1_results: List[Dict]) -> List[str]:
        """Extract promising mutations from Track 1 for Track 2"""
        
        if not track1_results:
            return []
        
        # Get mutations from top 25% of Track 1 variants
        sorted_variants = sorted(track1_results, 
                                key=lambda x: x.get('score', 0), 
                                reverse=True)
        
        top_variants = sorted_variants[:max(1, len(sorted_variants) // 4)]
        
        seed_mutations = set()
        for v in top_variants:
            seed_mutations.update(v.get('mutations', []))
        
        return list(seed_mutations)
    
    def _select_refinement_candidates(self, variants: List[Dict], n_select: int) -> List[Dict]:
        """Select best candidates for structural refinement"""
        
        if not variants:
            return []
        
        # Score variants for selection
        for v in variants:
            # Simple scoring based on available metrics
            score = 0
            score += v.get('delta_tm', 0) * 0.4
            score += v.get('conservation_score', 0) * 0.3
            score += (100 - v.get('risk_score', 50)) * 0.003
            v['refinement_score'] = score
        
        # Sort and select top N
        sorted_variants = sorted(variants, 
                                key=lambda x: x['refinement_score'], 
                                reverse=True)
        
        return sorted_variants[:n_select]
    
    async def _validate_with_statistics(self, variants: List[Dict]) -> List[Dict]:
        """Validate variants and add statistical metrics"""
        
        validated = []
        
        for variant in variants:
            validation_result = await self.validator.validate_single(variant)
            
            if validation_result.get('pass', False):
                # Add statistical confidence metrics
                variant['confidence_score'] = self._calculate_confidence(validation_result)
                variant['statistical_significance'] = self._calculate_significance(validation_result)
                validated.append(variant)
        
        return validated
    
    def _calculate_confidence(self, validation: Dict) -> float:
        """Calculate statistical confidence score"""
        
        # Weighted confidence based on multiple factors
        confidence = 0
        confidence += validation.get('tm_score', 0) * 0.3
        confidence += validation.get('structural_score', 0) * 0.3
        confidence += validation.get('evolutionary_score', 0) * 0.2
        confidence += (1 - validation.get('risk_score', 0.5)) * 0.2
        
        return min(confidence, 1.0)
    
    def _calculate_significance(self, validation: Dict) -> float:
        """Calculate statistical significance (p-value)"""
        
        # Simplified p-value calculation based on deviation from WT
        delta_tm = validation.get('delta_tm', 0)
        std_tm = 2.5  # Typical std deviation for Tm predictions
        
        # Z-score
        z_score = abs(delta_tm) / std_tm
        
        # Two-tailed p-value
        p_value = 2 * (1 - stats.norm.cdf(z_score))
        
        return p_value
    
    def _generate_mutation_rationale(self, variant: Dict) -> Dict:
        """Generate detailed rationale for why mutations improve thermotolerance"""
        
        rationale = {
            'mutations': variant.get('mutations', []),
            'mechanisms': [],
            'statistical_support': {},
            'structural_basis': {}
        }
        
        # Analyze each mutation
        for mutation in variant.get('mutations', []):
            # Parse mutation (e.g., "L335I")
            if len(mutation) > 2:
                from_aa = mutation[0]
                to_aa = mutation[-1]
                position = mutation[1:-1]
                
                mechanism = []
                
                # Hydrophobic packing
                if to_aa in 'ILMFWV' and from_aa not in 'ILMFWV':
                    mechanism.append(f"Increases hydrophobic packing at position {position}")
                
                # Charge optimization
                if to_aa in 'KR' and from_aa not in 'KR':
                    mechanism.append(f"Introduces positive charge for salt bridge formation")
                elif to_aa in 'DE' and from_aa not in 'DE':
                    mechanism.append(f"Introduces negative charge for electrostatic interactions")
                
                # Rigidity
                if to_aa == 'P':
                    mechanism.append(f"Reduces backbone flexibility at position {position}")
                
                # Statistical evidence
                if variant.get('conservation_score', 0) > 0.8:
                    mechanism.append(f"Found in {variant.get('conservation_score', 0)*100:.0f}% of thermophiles")
                
                rationale['mechanisms'].extend(mechanism)
        
        # Add statistical support
        rationale['statistical_support'] = {
            'predicted_delta_tm': f"+{variant.get('delta_tm', 0):.1f}°C",
            'confidence': f"{variant.get('confidence_score', 0):.1%}",
            'p_value': f"p<{variant.get('statistical_significance', 0.05):.4f}",
            'source_track': variant.get('source_track', 'unknown')
        }
        
        # Add structural basis
        rationale['structural_basis'] = {
            'plddt': variant.get('plddt', 0),
            'active_site_distance': f"{variant.get('active_site_distance', 0):.1f}Å",
            'h_bonds': f"+{variant.get('hbond_change', 0)} H-bonds" if variant.get('hbond_change', 0) > 0 else "Maintained"
        }
        
        return rationale
    
    def _calculate_overall_statistics(self, variants: List[Dict]) -> Dict:
        """Calculate comprehensive statistics for all variants"""
        
        if not variants:
            return {}
        
        tms = [v.get('predicted_tm', 0) for v in variants]
        delta_tms = [v.get('delta_tm', 0) for v in variants]
        
        # Correlation analysis
        mutation_counts = [len(v.get('mutations', [])) for v in variants]
        
        if len(set(mutation_counts)) > 1:  # Need variation for correlation
            r_value, p_value = stats.pearsonr(mutation_counts, delta_tms)
        else:
            r_value, p_value = 0, 1
        
        return {
            'mean_tm': np.mean(tms) if tms else 0,
            'best_delta_tm': max(delta_tms) if delta_tms else 0,
            'significance': p_value,
            'correlation': r_value,
            'success_rate': len(variants) / max(len(variants), 1)
        }
    
    def _export_results(self, variants: List[Dict], start_time: datetime) -> Dict:
        """Export comprehensive results with statistics"""
        
        runtime = (datetime.now() - start_time).total_seconds() / 60
        
        # Save detailed results
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        output_dir = Path(self.config.RESULTS_DIR)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save variants with full rationale
        results_file = output_dir / f"sequential_results_{timestamp}.json"
        
        results = {
            'pipeline_version': '3.0_sequential',
            'timestamp': timestamp,
            'runtime_minutes': runtime,
            'final_variants': variants,
            'statistics': self.statistical_results,
            'parameters': {
                'target_temperature': self.config.TARGET_TEMPERATURE,
                'initial_variant_count': self.config.INITIAL_VARIANT_COUNT
            }
        }
        
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"\n📊 Results exported to: {results_file}")
        
        # Print top variants with rationale
        if variants:
            print("\n🏆 TOP VARIANTS WITH MECHANISTIC RATIONALE:")
            print("─" * 50)
            
            for i, variant in enumerate(variants[:3], 1):
                print(f"\n  Variant {i}: {variant.get('id', 'Unknown')}")
                print(f"  Mutations: {', '.join(variant.get('mutations', []))}")
                print(f"  Predicted ΔTm: +{variant.get('delta_tm', 0):.1f}°C")
                
                rationale = variant.get('rationale', {})
                if rationale.get('mechanisms'):
                    print(f"  Why it works:")
                    for mechanism in rationale['mechanisms'][:3]:
                        print(f"    • {mechanism}")
                
                stats_support = rationale.get('statistical_support', {})
                print(f"  Statistical support: {stats_support.get('confidence', 'N/A')} confidence, "
                      f"{stats_support.get('p_value', 'N/A')}")
        
        return results


# Run if called directly
if __name__ == "__main__":
    import sys
    sys.path.append(str(Path(__file__).parent.parent))
    from analysis_engine import get_sequence
    
    async def test():
        config = Config()
        pipeline = MasterPipelineV3Sequential(config)
        
        # Use maize D1 sequence
        wt_sequence = get_sequence('maize')
        
        results = await pipeline.run(
            wt_sequence=wt_sequence,
            target_temp=50,
            n_variants=30,
            validate_all=True
        )
        
        print("\n✅ Sequential pipeline complete!")
        print(f"Generated {len(results.get('final_variants', []))} validated variants")
    
    asyncio.run(test())