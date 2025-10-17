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
from typing import List, Dict, Tuple
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
    
    def __init__(self, maize_sequence: str):
        """
        Args:
            maize_sequence: Wild-type maize D1 sequence
        """
        self.maize_sequence = maize_sequence
        
        # Protected regions for functional risk assessment
        self.active_sites = {130, 147, 161, 170, 189, 215, 254, 255, 264, 271, 332, 333, 342, 344}
    
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
            scored_variants.append(variant)
        
        # Sort by final score
        scored_variants.sort(key=lambda x: x['final_score'], reverse=True)
        
        print(f"   ✅ Scoring complete!")
        
        # Calculate statistics
        scores = [v['final_score'] for v in scored_variants]
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
                'mean_score': float(np.mean(scores)),
                'std_score': float(np.std(scores)),
                'max_score': float(np.max(scores)),
                'min_score': float(np.min(scores)),
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
            track = variant['source_track']
            for mutation in variant['mutations']:
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
        
        # Component 1: Thermostability (mock for now, would use TemStaPro API)
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
        for mutation in variant['mutations']:
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
        Estimate thermostability improvement (mock)
        
        In real implementation, use TemStaPro API
        """
        # Mock: base it on number of mutations and evidence score
        base = variant.get('evidence_score', 0.5)
        n_mut = variant['n_mutations']
        
        # More mutations = potentially more stability (up to a point)
        mut_factor = min(n_mut * 0.2, 0.6)
        
        # Random variation
        noise = np.random.uniform(-0.1, 0.1)
        
        return min(base + mut_factor + noise, 1.0)
    
    def _assess_evolutionary_evidence(self, variant: Dict) -> float:
        """
        Assess evolutionary support for variant
        """
        track = variant['source_track']
        
        # Evolutionary track has strongest evidence
        if track == 'evolutionary_guided':
            return variant.get('evidence_score', 0.8)
        
        # Other tracks have moderate evidence
        return 0.5
    
    def _assess_functional_risk(self, variant: Dict) -> float:
        """
        Assess risk of disrupting protein function
        
        Risk = proximity to active sites
        """
        min_distance = 100.0
        
        for mutation in variant['mutations']:
            # Parse mutation (e.g., "L335I")
            pos = int(''.join(filter(str.isdigit, mutation)))
            
            # Distance to nearest active site
            distance = min(abs(pos - site) for site in self.active_sites)
            min_distance = min(min_distance, distance)
        
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
        
        if len(variants) <= n_clusters:
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
        pos1 = set(int(''.join(filter(str.isdigit, m))) for m in var1['mutations'])
        pos2 = set(int(''.join(filter(str.isdigit, m))) for m in var2['mutations'])
        
        # Symmetric difference
        return len(pos1.symmetric_difference(pos2))
    
    def export_results(self, output_dir: Path, results: Dict):
        """Export consensus results"""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Export top variants as FASTA
        fasta_file = output_dir / "final_top_variants.fasta"
        with open(fasta_file, 'w') as f:
            for i, variant in enumerate(results['top_variants'], 1):
                mutations_str = "+".join(variant['mutations'])
                f.write(f">{variant['variant_id']} | Rank={i} | {mutations_str} | "
                       f"Score={variant['final_score']:.3f} | "
                       f"Consensus={'YES' if variant['has_consensus'] else 'NO'}\n")
                seq = variant['sequence']
                for j in range(0, len(seq), 60):
                    f.write(f"{seq[j:j+60]}\n")
        
        print(f"\n💾 Top variants FASTA: {fasta_file}")
        
        # Export scoring table as CSV
        csv_file = output_dir / "variant_scores.csv"
        df_data = []
        for variant in results['scored_variants']:
            df_data.append({
                'variant_id': variant['variant_id'],
                'mutations': '+'.join(variant['mutations']),
                'n_mutations': variant['n_mutations'],
                'source_track': variant['source_track'],
                'final_score': variant['final_score'],
                'thermostability': variant['thermostability_component'],
                'structure': variant['structure_component'],
                'evolutionary': variant['evolutionary_component'],
                'functional_risk': variant['functional_risk_component'],
                'has_consensus': variant['has_consensus'],
                'identity_to_wt': variant['identity_to_wt']
            })
        
        df = pd.DataFrame(df_data)
        df.to_csv(csv_file, index=False)
        print(f"💾 Scoring table: {csv_file}")
        
        # Export consensus mutations
        consensus_file = output_dir / "consensus_mutations.json"
        with open(consensus_file, 'w') as f:
            json.dump(results['consensus_mutations'], f, indent=2)
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
    from analysis_engine import REFERENCE_SEQUENCES
    
    print("="*70)
    print("🏆 CONSENSUS SCORER TEST")
    print("="*70)
    print("Loading mock results from all tracks...")
    
    # Mock results for testing
    maize_seq = REFERENCE_SEQUENCES['maize']
    
    track1_results = {'variants': []}  # Track 1 had 0 variants
    track2_results = {'variants': [
        {
            'variant_id': 'TRACK2_EE_001',
            'sequence': maize_seq,
            'mutations': ['L335I'],
            'n_mutations': 1,
            'identity_to_wt': 0.997,
            'evidence_score': 0.75,
            'source_track': 'generative_ai'
        }
    ]}
    track3_results = {'variants': [
        {
            'variant_id': 'TRACK3_STRUCT_001',
            'sequence': maize_seq,
            'mutations': ['L335I'],
            'n_mutations': 1,
            'identity_to_wt': 0.997,
            'evidence_score': 0.82,
            'source_track': 'structure_guided',
            'structural_confidence': 0.89
        }
    ]}
    
    # Run consensus scoring
    scorer = ConsensusScorer(maize_seq)
    results = scorer.score_all_variants(track1_results, track2_results, track3_results)
    
    # Export
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(f"../data/results/consensus_{timestamp}")
    scorer.export_results(output_dir, results)
    
    print("\n" + "="*70)
    print("✅ CONSENSUS SCORING COMPLETE!")
    print("="*70)
    print(f"📊 Scored {results['total_variants']} variants")
    print(f"🏆 Consensus mutations: {len(results['consensus_mutations'])}")
    if results['top_variants']:
        print(f"🥇 Top variant: {results['top_variants'][0]['variant_id']} "
              f"(Score: {results['top_variants'][0]['final_score']:.2f})")
    print(f"📁 Results saved to: {output_dir}")
    print("="*70)