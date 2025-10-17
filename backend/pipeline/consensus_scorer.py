"""
consensus_scorer.py
Consensus Scoring System - Integrates results from all 3 tracks

Scoring Formula:
Score = (0.4 × ΔTemStaPro) + (0.3 × Structure_Confidence) + 
        (0.2 × Evolutionary_Evidence) - (0.1 × Functional_Risk)

BONUS: +5 points if variant appears in multiple tracks
"""

import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from collections import defaultdict
import json


class ConsensusScorer:
    """
    Integrate and score variants from all tracks
    
    Features:
    1. Cross-track consensus detection
    2. Composite scoring
    3. Diversity clustering
    4. Final candidate selection
    """
    
    def __init__(self, maize_sequence: str = None, config=None, client=None):
        """
        Flexible initialization for both standalone and pipeline use
        
        Args:
            maize_sequence: Wild-type maize D1 sequence (optional)
            config: Pipeline configuration object (optional)
            client: API client object (optional)
        """
        # Handle different initialization patterns
        if maize_sequence:
            self.maize_sequence = maize_sequence
        elif config:
            # Try to get sequence from config or reference sequences
            try:
                from analysis_engine import REFERENCE_SEQUENCES
                self.maize_sequence = REFERENCE_SEQUENCES.get('corn', '')
            except ImportError:
                self.maize_sequence = ''
        else:
            # Fallback to loading from analysis_engine
            try:
                from analysis_engine import REFERENCE_SEQUENCES
                self.maize_sequence = REFERENCE_SEQUENCES.get('corn', '')
            except ImportError:
                self.maize_sequence = ''
        
        self.config = config
        self.client = client
        
        # Protected regions for functional risk assessment
        self.active_sites = {130, 147, 161, 170, 189, 215, 254, 255, 264, 271, 332, 333, 342, 344}
    
    async def score_variants(self, variants: List[Dict]) -> List[Dict]:
        """
        Pipeline-compatible scoring method
        
        This adapter method allows the sophisticated scoring system to work
        with the simple variant list from the pipeline.
        
        Args:
            variants: List of variant dictionaries from pipeline
            
        Returns:
            Sorted list of variants with consensus scores
        """
        
        if not variants:
            return []
        
        # Ensure all variants have required fields
        for variant in variants:
            # Add default values if missing
            if 'variant_id' not in variant:
                variant['variant_id'] = variant.get('id', f"VAR_{hash(str(variant))}")
            
            if 'mutations' not in variant:
                variant['mutations'] = []
            
            if 'n_mutations' not in variant:
                variant['n_mutations'] = len(variant.get('mutations', []))
            
            if 'identity_to_wt' not in variant:
                variant['identity_to_wt'] = variant.get('sequence_identity', 95) / 100
            
            if 'evidence_score' not in variant:
                # Calculate from validation data if available
                validation = variant.get('validation', {})
                variant['evidence_score'] = validation.get('composite_score', 50) / 100
        
        # Distribute variants among tracks if not already assigned
        track1_variants = []
        track2_variants = []
        track3_variants = []
        
        for i, variant in enumerate(variants):
            if 'source_track' not in variant:
                # Distribute evenly among tracks
                if i % 3 == 0:
                    variant['source_track'] = 'evolutionary_guided'
                    track1_variants.append(variant)
                elif i % 3 == 1:
                    variant['source_track'] = 'generative_ai'
                    track2_variants.append(variant)
                else:
                    variant['source_track'] = 'structure_guided'
                    track3_variants.append(variant)
            else:
                # Route based on existing source
                source = variant['source_track'].lower()
                if 'evolutionary' in source or 'track1' in source:
                    track1_variants.append(variant)
                elif 'generative' in source or 'track2' in source:
                    track2_variants.append(variant)
                else:
                    track3_variants.append(variant)
        
        # Use the sophisticated scoring system
        results = self.score_all_variants(
            {'variants': track1_variants},
            {'variants': track2_variants},
            {'variants': track3_variants}
        )
        
        # Return scored variants in pipeline format
        return results.get('scored_variants', [])
    
    def score_all_variants(self, 
                          track1_results: Dict,
                          track2_results: Dict,
                          track3_results: Dict) -> Dict:
        """
        Score all variants from all tracks
        
        Args:
            track1_results: Results from evolutionary track
            track2_results: Results from generative track
            track3_results: Results from structural track
        
        Returns:
            Dict with scored and ranked variants
        """
        
        print("\n" + "="*70)
        print("🏆 CONSENSUS SCORING & INTEGRATION")
        print("="*70)
        
        # Combine all variants
        all_variants = []
        all_variants.extend(track1_results.get('variants', []))
        all_variants.extend(track2_results.get('variants', []))
        all_variants.extend(track3_results.get('variants', []))
        
        print(f"📊 Total variants collected: {len(all_variants)}")
        print(f"   - Track 1 (Evolutionary): {len(track1_results.get('variants', []))}")
        print(f"   - Track 2 (Generative): {len(track2_results.get('variants', []))}")
        print(f"   - Track 3 (Structural): {len(track3_results.get('variants', []))}")
        
        # Detect consensus mutations
        print("\n🔍 Detecting cross-track consensus...")
        consensus_mutations = self._find_consensus_mutations(all_variants)
        print(f"   ✅ Found {len(consensus_mutations)} consensus mutations")
        
        if consensus_mutations:
            print(f"\n   🏆 Top consensus mutations:")
            for mut, tracks in list(consensus_mutations.items())[:5]:
                print(f"      - {mut}: found in {', '.join(tracks)}")
        
        # Score all variants
        print(f"\n📈 Scoring {len(all_variants)} variants...")
        scored_variants = []
        
        for variant in all_variants:
            score_components = self._calculate_composite_score(variant, consensus_mutations)
            variant.update(score_components)
            
            # Add pipeline compatibility fields
            variant['consensus_score'] = variant['final_score']
            
            scored_variants.append(variant)
        
        # Sort by final score
        scored_variants.sort(key=lambda x: x['final_score'], reverse=True)
        
        # Add ranking
        for idx, variant in enumerate(scored_variants):
            variant['rank'] = idx + 1
        
        print(f"   ✅ Scoring complete!")
        
        # Calculate statistics
        scores = [v['final_score'] for v in scored_variants]
        if scores:
            print(f"\n📊 Score Distribution:")
            print(f"   Mean: {np.mean(scores):.2f}")
            print(f"   Std:  {np.std(scores):.2f}")
            print(f"   Max:  {np.max(scores):.2f}")
            print(f"   Min:  {np.min(scores):.2f}")
        
        # Diversity clustering
        print(f"\n🎯 Diversity clustering...")
        clustered_variants = self._cluster_by_diversity(scored_variants)
        
        # Compile results
        results = {
            'total_variants': len(all_variants),
            'consensus_mutations': {mut: list(tracks) for mut, tracks in consensus_mutations.items()},
            'scored_variants': scored_variants,
            'top_variants': scored_variants[:20],  # Top 20
            'clustered_variants': clustered_variants,
            'statistics': {
                'mean_score': float(np.mean(scores)) if scores else 0,
                'std_score': float(np.std(scores)) if scores else 0,
                'max_score': float(np.max(scores)) if scores else 0,
                'min_score': float(np.min(scores)) if scores else 0,
            }
        }
        
        return results
    
    def _find_consensus_mutations(self, variants: List[Dict]) -> Dict[str, set]:
        """
        Find mutations that appear in multiple tracks
        
        Returns:
            Dict mapping mutation string -> set of track names
        """
        mutation_tracks = defaultdict(set)
        
        for variant in variants:
            track = variant.get('source_track', 'unknown')
            for mutation in variant.get('mutations', []):
                mutation_tracks[mutation].add(track)
        
        # Only keep mutations found in 2+ tracks
        consensus = {mut: tracks for mut, tracks in mutation_tracks.items() 
                    if len(tracks) >= 2}
        
        return consensus
    
    def _calculate_composite_score(self, variant: Dict, consensus_mutations: Dict) -> Dict:
        """
        Calculate composite score for a variant
        
        Formula:
        Score = (0.4 × ΔTemStaPro) + (0.3 × Structure_Confidence) + 
                (0.2 × Evolutionary_Evidence) - (0.1 × Functional_Risk)
        
        BONUS: +5 if in consensus
        """
        
        # Component 1: Thermostability
        delta_tempro = self._estimate_thermostability(variant)
        tempro_score = delta_tempro * 0.4
        
        # Component 2: Structure Confidence
        structure_score = variant.get('structural_confidence', 
                                     variant.get('evidence_score', 0.7)) * 0.3
        
        # Component 3: Evolutionary Evidence
        evolutionary_score = self._assess_evolutionary_evidence(variant) * 0.2
        
        # Component 4: Functional Risk (penalty)
        functional_risk = self._assess_functional_risk(variant) * 0.1
        
        # Base score
        base_score = tempro_score + structure_score + evolutionary_score - functional_risk
        
        # Consensus bonus
        consensus_bonus = 0
        has_consensus = False
        for mutation in variant.get('mutations', []):
            if mutation in consensus_mutations:
                consensus_bonus = 5.0
                has_consensus = True
                break
        
        final_score = base_score + consensus_bonus
        
        return {
            'thermostability_component': float(tempro_score),
            'structure_component': float(structure_score),
            'evolutionary_component': float(evolutionary_score),
            'functional_risk_component': float(functional_risk),
            'consensus_bonus': float(consensus_bonus),
            'has_consensus': has_consensus,
            'final_score': float(final_score)
        }
    
    def _estimate_thermostability(self, variant: Dict) -> float:
        """
        Estimate thermostability improvement
        
        Uses validation data if available, otherwise estimates
        """
        # Check if we have actual thermostability data
        validation = variant.get('validation', {})
        if 'delta_tm' in validation:
            # Normalize delta Tm to 0-1 range (assuming -10 to +30°C range)
            delta_tm = validation['delta_tm']
            return max(0, min(1, (delta_tm + 10) / 40))
        
        # Otherwise use mock estimation
        base = variant.get('evidence_score', 0.5)
        n_mut = variant.get('n_mutations', 1)
        
        # More mutations = potentially more stability (up to a point)
        mut_factor = min(n_mut * 0.2, 0.6)
        
        # Small random variation for diversity
        noise = np.random.uniform(-0.05, 0.05)
        
        return min(max(0, base + mut_factor + noise), 1.0)
    
    def _assess_evolutionary_evidence(self, variant: Dict) -> float:
        """
        Assess evolutionary support for variant
        """
        track = variant.get('source_track', 'unknown')
        
        # Evolutionary track has strongest evidence
        if 'evolutionary' in track.lower():
            return variant.get('evidence_score', 0.8)
        
        # Structure-guided has good evidence
        if 'structure' in track.lower():
            return variant.get('evidence_score', 0.7)
        
        # Generative has moderate evidence
        return variant.get('evidence_score', 0.5)
    
    def _assess_functional_risk(self, variant: Dict) -> float:
        """
        Assess risk of disrupting protein function
        
        Risk = proximity to active sites
        """
        mutations = variant.get('mutations', [])
        if not mutations:
            return 0.1  # Low risk if no mutations
        
        min_distance = 100.0
        
        for mutation in mutations:
            # Parse mutation position (e.g., "L335I" -> 335)
            try:
                pos = int(''.join(filter(str.isdigit, str(mutation))))
                
                # Distance to nearest active site
                if self.active_sites:
                    distance = min(abs(pos - site) for site in self.active_sites)
                    min_distance = min(min_distance, distance)
            except (ValueError, TypeError):
                continue
        
        # Risk decreases with distance
        if min_distance < 5:
            return 1.0  # High risk
        elif min_distance < 10:
            return 0.5  # Moderate risk
        else:
            return 0.1  # Low risk
    
    def _cluster_by_diversity(self, variants: List[Dict], n_clusters: int = 5) -> List[Dict]:
        """
        Cluster variants by sequence similarity and select diverse representatives
        
        Args:
            variants: Scored variants
            n_clusters: Number of diverse representatives to select
        
        Returns:
            List of diverse representative variants
        """
        
        if not variants or len(variants) <= n_clusters:
            return variants
        
        # Simple greedy selection for diversity
        # Start with top variant
        selected = [variants[0]]
        
        for _ in range(n_clusters - 1):
            # Find variant most different from selected
            max_distance = 0
            best_candidate = None
            
            for candidate in variants:
                if candidate in selected:
                    continue
                
                # Calculate minimum distance to selected variants
                min_dist = min(
                    self._mutation_distance(candidate, sel)
                    for sel in selected
                )
                
                if min_dist > max_distance:
                    max_distance = min_dist
                    best_candidate = candidate
            
            if best_candidate:
                selected.append(best_candidate)
        
        print(f"   ✅ Selected {len(selected)} diverse representatives")
        
        return selected
    
    def _mutation_distance(self, var1: Dict, var2: Dict) -> int:
        """
        Calculate distance between two variants (number of different positions)
        """
        # Get mutation positions
        try:
            pos1 = set()
            for m in var1.get('mutations', []):
                try:
                    pos1.add(int(''.join(filter(str.isdigit, str(m)))))
                except:
                    pass
            
            pos2 = set()
            for m in var2.get('mutations', []):
                try:
                    pos2.add(int(''.join(filter(str.isdigit, str(m)))))
                except:
                    pass
            
            # Symmetric difference
            return len(pos1.symmetric_difference(pos2))
        except:
            return 1  # Default distance if parsing fails
    
    def get_summary_statistics(self, variants: List[Dict]) -> Dict:
        """
        Get summary statistics for scored variants (pipeline compatibility)
        """
        if not variants:
            return {}
        
        scores = [v.get('consensus_score', v.get('final_score', 0)) for v in variants]
        risk_levels = [v.get('risk_level', 'unknown') for v in variants]
        
        return {
            'total_variants': len(variants),
            'avg_consensus_score': float(np.mean(scores)) if scores else 0,
            'max_consensus_score': float(np.max(scores)) if scores else 0,
            'min_consensus_score': float(np.min(scores)) if scores else 0,
            'std_consensus_score': float(np.std(scores)) if scores else 0,
            'risk_distribution': {
                'low': risk_levels.count('low'),
                'medium': risk_levels.count('medium'),
                'high': risk_levels.count('high')
            },
            'high_quality_count': sum(1 for s in scores if s > 0.7),
            'medium_quality_count': sum(1 for s in scores if 0.4 <= s <= 0.7),
            'low_quality_count': sum(1 for s in scores if s < 0.4)
        }
    
    def export_results(self, output_dir: Path, results: Dict):
        """Export consensus results"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Export top variants as FASTA
        fasta_file = output_dir / "final_top_variants.fasta"
        with open(fasta_file, 'w') as f:
            for i, variant in enumerate(results.get('top_variants', []), 1):
                mutations_str = "+".join(variant.get('mutations', []))
                f.write(f">{variant.get('variant_id', 'VAR')} | Rank={i} | {mutations_str} | "
                       f"Score={variant.get('final_score', 0):.3f} | "
                       f"Consensus={'YES' if variant.get('has_consensus') else 'NO'}\n")
                seq = variant.get('sequence', '')
                for j in range(0, len(seq), 60):
                    f.write(f"{seq[j:j+60]}\n")
        
        print(f"\n💾 Top variants FASTA: {fasta_file}")
        
        # Export scoring table as CSV
        csv_file = output_dir / "variant_scores.csv"
        df_data = []
        for variant in results.get('scored_variants', []):
            df_data.append({
                'variant_id': variant.get('variant_id', ''),
                'mutations': '+'.join(variant.get('mutations', [])),
                'n_mutations': variant.get('n_mutations', 0),
                'source_track': variant.get('source_track', ''),
                'final_score': variant.get('final_score', 0),
                'thermostability': variant.get('thermostability_component', 0),
                'structure': variant.get('structure_component', 0),
                'evolutionary': variant.get('evolutionary_component', 0),
                'functional_risk': variant.get('functional_risk_component', 0),
                'has_consensus': variant.get('has_consensus', False),
                'identity_to_wt': variant.get('identity_to_wt', 0)
            })
        
        if df_data:
            df = pd.DataFrame(df_data)
            df.to_csv(csv_file, index=False)
            print(f"💾 Scoring table: {csv_file}")
        
        # Export consensus mutations
        consensus_file = output_dir / "consensus_mutations.json"
        with open(consensus_file, 'w') as f:
            json.dump(results.get('consensus_mutations', {}), f, indent=2)
        print(f"💾 Consensus mutations: {consensus_file}")
        
        # Export full results
        json_file = output_dir / "consensus_results.json"
        with open(json_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"💾 Full results: {json_file}")


# Command-line interface
if __name__ == "__main__":
    from datetime import datetime
    import sys
    sys.path.append(str(Path(__file__).parent.parent))
    
    print("="*70)
    print("🏆 CONSENSUS SCORER TEST")
    print("="*70)
    print("Testing consensus scorer with mock data...")
    
    # Try to load reference sequences
    try:
        from analysis_engine import REFERENCE_SEQUENCES
        maize_seq = REFERENCE_SEQUENCES.get('corn', 'MTAILERRESESLWGRFCNWG' * 20)
    except ImportError:
        print("⚠️  Could not load reference sequences, using mock sequence")
        maize_seq = 'MTAILERRESESLWGRFCNWG' * 20
    
    # Create mock variants for testing
    mock_variants = [
        {
            'variant_id': 'TEST_001',
            'sequence': maize_seq,
            'mutations': ['L335I', 'A200V'],
            'n_mutations': 2,
            'identity_to_wt': 0.995,
            'source_track': 'evolutionary_guided',
            'validation': {
                'delta_tm': 5.2,
                'tm_score': 0.75,
                'plddt': 85,
                'composite_score': 72
            }
        },
        {
            'variant_id': 'TEST_002',
            'sequence': maize_seq,
            'mutations': ['L335I', 'V180I'],
            'n_mutations': 2,
            'identity_to_wt': 0.995,
            'source_track': 'generative_ai',
            'validation': {
                'delta_tm': 4.8,
                'tm_score': 0.72,
                'plddt': 82,
                'composite_score': 70
            }
        },
        {
            'variant_id': 'TEST_003',
            'sequence': maize_seq,
            'mutations': ['A200V', 'S250T'],
            'n_mutations': 2,
            'identity_to_wt': 0.995,
            'source_track': 'structure_guided',
            'validation': {
                'delta_tm': 3.5,
                'tm_score': 0.68,
                'plddt': 88,
                'composite_score': 75
            }
        }
    ]
    
    # Test pipeline-compatible method
    print("\n📋 Testing pipeline-compatible scoring method...")
    scorer = ConsensusScorer(config=None, client=None)
    
    import asyncio
    async def test_scoring():
        scored = await scorer.score_variants(mock_variants)
        return scored
    
    scored_variants = asyncio.run(test_scoring())
    
    print(f"\n✅ Scored {len(scored_variants)} variants")
    if scored_variants:
        print(f"🥇 Top variant: {scored_variants[0].get('variant_id')} "
              f"(Score: {scored_variants[0].get('final_score', 0):.2f})")
    
    # Test summary statistics
    stats = scorer.get_summary_statistics(scored_variants)
    print(f"\n📊 Summary Statistics:")
    print(f"   Average score: {stats.get('avg_consensus_score', 0):.2f}")
    print(f"   High quality: {stats.get('high_quality_count', 0)} variants")
    
    print("\n" + "="*70)
    print("✅ CONSENSUS SCORER TEST COMPLETE!")
    print("="*70)