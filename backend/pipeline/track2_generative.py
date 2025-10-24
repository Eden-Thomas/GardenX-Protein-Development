"""
track2_generative.py
AI-guided optimization with multi-objective scoring and Track 1 seed integration
UPDATED: Added iteratively_improve_variant() for true sequential pipeline
"""

import asyncio
import numpy as np
from typing import Dict, List, Optional
import random


class Track2Generative:
    """AI-guided variant generation building on Track 1 insights"""
    
    def __init__(self, config, client):
        self.config = config
        self.client = client
        self.membrane_regions = {
            'TM1': range(21, 47),
            'TM2': range(60, 88),
            'TM3': range(181, 207),
            'TM4': range(241, 270),
            'TM5': range(311, 342)
        }
        
        self.functional_exclusion = {170, 189, 161, 332, 333, 342, 344, 
                                    215, 254, 255, 264, 265, 271, 130, 147}
    
    async def generate_variants(
        self,
        wt_sequence: str,
        target_temp: int = 50,
        n_variants: int = 15,
        seed_mutations: List[str] = None,
        seed_variants: List[Dict] = None
    ) -> List[Dict]:
        """
        Generate AI-optimized variants using Track 1 insights
        
        Args:
            wt_sequence: Wild-type sequence
            target_temp: Target temperature (°C)
            n_variants: Number of variants to generate
            seed_mutations: Top mutations from Track 1
            seed_variants: Full variants from Track 1 (for reference)
        """
        
        print("  🤖 Track 2: AI-guided optimization")
        
        if seed_mutations and len(seed_mutations) > 0:
            print(f"    → Building on {len(seed_mutations)} Track 1 mutations")
            print(f"      Seeds: {', '.join(seed_mutations[:5])}...")
        else:
            print("    → No Track 1 seeds, generating de novo")
            seed_mutations = []
        
        variants = []
        
        # Strategy 1: NeuroFold-based design (if API available)
        if self.client and len(variants) < n_variants * 0.3:
            try:
                print("    → Trying NeuroFold design...")
                neurofold_variants = await self._neurofold_with_constraints(
                    wt_sequence,
                    seed_mutations,
                    max(1, int(n_variants * 0.3))
                )
                if neurofold_variants:
                    variants.extend(neurofold_variants)
                    print(f"      ✓ NeuroFold generated {len(neurofold_variants)} variants")
            except Exception as e:
                print(f"      ⚠️ NeuroFold failed: {e}")
        
        # Strategy 2: Track 1-guided fallback (PRIMARY METHOD)
        if seed_mutations and len(variants) < n_variants:
            print("    → Applying Track 1-guided mutations...")
            remaining = n_variants - len(variants)
            track1_variants = self._apply_track1_mutations(
                wt_sequence,
                seed_mutations,
                remaining
            )
            variants.extend(track1_variants)
            print(f"      ✓ Generated {len(track1_variants)} Track 1-guided variants")
        
        # Strategy 3: Pure fallback if still short
        if len(variants) < n_variants:
            print("    → Filling remaining with heuristic variants...")
            remaining = n_variants - len(variants)
            fallback = self._heuristic_fallback(wt_sequence, remaining)
            variants.extend(fallback)
        
        # Multi-objective scoring for all variants
        print("    → Scoring variants (multi-objective)...")
        for i, variant in enumerate(variants):
            scores = await self._multi_objective_score(variant, wt_sequence, target_temp)
            variant.update(scores)
            if 'id' not in variant:
                variant['id'] = f'T2_AI_{i:04d}'
            variant['track'] = 'track2'
        
        variants.sort(key=lambda x: x.get('composite_score', 0), reverse=True)
        
        print(f"    ✓ Generated {len(variants)} AI-optimized variants")
        if variants:
            print(f"      Mean ΔTm: {np.mean([v.get('delta_tm', 0) for v in variants]):.1f}°C")
        
        return variants[:n_variants]
    
    async def iteratively_improve_variant(
        self,
        base_variant: Dict,
        wt_sequence: str,
        target_temp: int = 55,
        improvement_rounds: int = 3,
        aggressive_mode: bool = True
    ) -> Dict:
        """
        NEW METHOD: Take an EXACT variant from Track 1 and iteratively improve it
        This is the key method for true sequential processing
        
        Args:
            base_variant: Variant from Track 1 to improve
            wt_sequence: Wild-type sequence
            target_temp: Target temperature
            improvement_rounds: Number of iterative improvement rounds
            aggressive_mode: Use aggressive thermostabilizing mutations
        """
        
        print(f"      Iteratively improving {base_variant['id']} (target: {target_temp}°C)")
        
        current_variant = base_variant.copy()
        current_tm = current_variant.get('predicted_tm', 28)
        
        for round_num in range(improvement_rounds):
            print(f"        Round {round_num + 1}: Current Tm = {current_tm:.1f}°C")
            
            # Find positions near existing mutations for synergistic effects
            existing_positions = set()
            for mut in current_variant.get('mutations', []):
                if len(mut) >= 3:
                    try:
                        pos = int(mut[1:-1])
                        existing_positions.add(pos)
                    except:
                        continue
            
            # Generate candidate positions (within 3 residues of existing mutations)
            candidate_positions = set()
            for pos in existing_positions:
                for offset in range(-3, 4):  # ±3 residues
                    new_pos = pos + offset
                    if 10 <= new_pos <= len(wt_sequence) - 10:  # Safe range
                        candidate_positions.add(new_pos)
            
            # Add some random positions for diversity
            for _ in range(5):
                random_pos = np.random.randint(50, len(wt_sequence) - 50)
                candidate_positions.add(random_pos)
            
            # Generate thermostabilizing mutations at these positions
            improvement_mutations = []
            for pos in list(candidate_positions)[:10]:  # Limit to 10 candidates
                if pos in self.functional_exclusion:
                    continue
                
                pos_idx = pos - 1
                if pos_idx >= len(wt_sequence):
                    continue
                
                wt_aa = wt_sequence[pos_idx]
                
                # Aggressive thermostabilizing substitutions
                if aggressive_mode:
                    thermostable_subs = {
                        'G': ['A', 'V', 'I'],      # Remove flexibility
                        'S': ['T', 'V', 'A'],      # Reduce polarity
                        'N': ['D', 'Q', 'T'],      # Improve packing
                        'Q': ['E', 'N', 'L'],      # Enhance stability
                        'A': ['V', 'I', 'L'],      # Increase hydrophobicity
                        'T': ['V', 'I', 'A'],      # Better packing
                        'V': ['I', 'L'],           # Increase size
                        'L': ['I', 'V'],           # Optimize packing
                    }
                else:
                    # Conservative substitutions
                    thermostable_subs = {
                        'G': ['A'],
                        'S': ['T'],
                        'A': ['V']
                    }
                
                if wt_aa in thermostable_subs:
                    to_aa = np.random.choice(thermostable_subs[wt_aa])
                    improvement_mutations.append(f'{wt_aa}{pos}{to_aa}')
            
            # Select 2-4 improvement mutations
            if improvement_mutations:
                n_improvements = np.random.randint(2, min(5, len(improvement_mutations) + 1))
                selected_improvements = np.random.choice(
                    improvement_mutations, 
                    size=n_improvements, 
                    replace=False
                )
                
                # Apply improvements
                new_mutations = current_variant['mutations'] + list(selected_improvements)
                new_sequence = self._apply_mutations_to_sequence(wt_sequence, new_mutations)
                
                # Predict new Tm (use aggressive prediction)
                if aggressive_mode:
                    new_tm = self._aggressive_tm_prediction(new_sequence, wt_sequence)
                else:
                    new_tm = current_tm + len(selected_improvements) * 1.5  # Heuristic
                
                # Keep if improvement
                if new_tm > current_tm:
                    current_variant = {
                        'id': f"{base_variant['id']}_R{round_num + 1}",
                        'sequence': new_sequence,
                        'mutations': new_mutations,
                        'predicted_tm': new_tm,
                        'track': 'track2',
                        'method': 'iterative_improvement',
                        'track1_parent': base_variant['id'],
                        'improvement_round': round_num + 1,
                        'improvements_added': list(selected_improvements)
                    }
                    current_tm = new_tm
                    print(f"          Improved to {new_tm:.1f}°C (+{len(selected_improvements)} mutations)")
                else:
                    print(f"          No improvement in round {round_num + 1}")
            else:
                print(f"          No candidate improvements found")
        
        # Final variant
        final_variant = current_variant.copy()
        final_variant['id'] = f"{base_variant['id']}_IMPROVED"
        
        print(f"      Final: {base_variant['id']} → {final_variant['id']}")
        print(f"      Tm: {base_variant.get('predicted_tm', 28):.1f}°C → {final_variant.get('predicted_tm', 28):.1f}°C")
        print(f"      Mutations: {len(base_variant.get('mutations', []))} → {len(final_variant.get('mutations', []))}")
        
        return final_variant
    
    def _aggressive_tm_prediction(self, sequence: str, wt_sequence: str) -> float:
        """
        NEW METHOD: More aggressive Tm prediction targeting higher temperatures
        """
        
        base_tm = 28.0  # WT baseline
        
        # Count different types of stabilizing changes
        mutations = self._get_mutations(wt_sequence, sequence)
        
        tm_boost = 0.0
        for mut in mutations:
            if len(mut) < 3:
                continue
            
            from_aa = mut[0]
            to_aa = mut[-1]
            
            # Higher scoring for known thermostable substitutions
            if from_aa == 'G' and to_aa in 'AVIL':
                tm_boost += 3.0  # Large boost for Gly→hydrophobic
            elif from_aa == 'S' and to_aa in 'TVA':
                tm_boost += 2.5  # Good boost for Ser→less polar
            elif from_aa in 'NQ' and to_aa in 'DETL':
                tm_boost += 2.0  # Moderate boost for polar optimization
            elif from_aa in 'AV' and to_aa in 'IL':
                tm_boost += 1.5  # Small boost for size increase
            else:
                tm_boost += 1.0  # Generic boost
        
        # Diminishing returns
        if tm_boost > 25:
            tm_boost = 25 + (tm_boost - 25) * 0.5
        
        return base_tm + tm_boost
    
    async def _neurofold_with_constraints(
        self, wt_sequence: str, seed_mutations: List[str], n_variants: int
    ) -> List[Dict]:
        """Use NeuroFold with Track 1 mutations as hints"""
        
        variants = []
        
        if seed_mutations:
            template = self._apply_mutations_to_sequence(wt_sequence, seed_mutations[:3])
        else:
            template = wt_sequence
        
        try:
            result = await self.client.run_neurofold(template)
            
            if 'sequences' in result and result['sequences']:
                for i, seq_data in enumerate(result['sequences'][:n_variants]):
                    sequence = seq_data.get('sequence', seq_data) if isinstance(seq_data, dict) else seq_data
                    mutations = self._get_mutations(wt_sequence, sequence)
                    
                    variants.append({
                        'id': f'T2_NEUROFOLD_{i:04d}',
                        'sequence': sequence,
                        'mutations': mutations,
                        'method': 'neurofold',
                        'based_on_track1': len(seed_mutations) > 0
                    })
        except Exception as e:
            print(f"      ⚠️ NeuroFold error: {e}")
        
        return variants
    
    def _apply_track1_mutations(
        self, wt_sequence: str, seed_mutations: List[str], n_variants: int
    ) -> List[Dict]:
        """
        Generate variants by combining Track 1 mutations with synergistic additions
        THIS IS THE KEY METHOD FOR SEQUENTIAL PIPELINE!
        """
        
        variants = []
        
        parsed_seeds = []
        for mut in seed_mutations[:10]:
            if len(mut) >= 3:
                try:
                    from_aa = mut[0]
                    pos = int(mut[1:-1])
                    to_aa = mut[-1]
                    parsed_seeds.append({'from': from_aa, 'pos': pos, 'to': to_aa, 'string': mut})
                except:
                    continue
        
        if not parsed_seeds:
            return []
        
        for i in range(n_variants):
            n_t1_muts = min(np.random.randint(3, 6), len(parsed_seeds))
            selected_t1 = random.sample(parsed_seeds, n_t1_muts)
            
            variant_seq = list(wt_sequence)
            variant_muts = []
            
            for mut in selected_t1:
                pos_idx = mut['pos'] - 1
                if pos_idx < len(variant_seq) and variant_seq[pos_idx] == mut['from']:
                    variant_seq[pos_idx] = mut['to']
                    variant_muts.append(mut['string'])
            
            synergistic_muts = self._add_synergistic_mutations(
                variant_seq, wt_sequence, selected_t1, np.random.randint(1, 4)
            )
            
            variant_muts.extend(synergistic_muts)
            
            variants.append({
                'id': f'T2_FALLBACK_{i:04d}',
                'sequence': ''.join(variant_seq),
                'mutations': variant_muts,
                'method': 'track1_guided',
                'n_track1_muts': len(selected_t1),
                'n_synergistic': len(synergistic_muts),
                'based_on_track1': True
            })
        
        return variants
    
    def _add_synergistic_mutations(
        self, variant_seq: List[str], wt_sequence: str, 
        existing_muts: List[Dict], n_synergistic: int = 2
    ) -> List[str]:
        """Add synergistic mutations near Track 1 sites"""
        
        synergistic = []
        existing_positions = {m['pos'] for m in existing_muts}
        
        stabilizing = {
            'G': ['A', 'V'],
            'S': ['T', 'A'],
            'N': ['D', 'Q'],
            'Q': ['E', 'N']
        }
        
        attempts = 0
        while len(synergistic) < n_synergistic and attempts < 50:
            attempts += 1
            
            if existing_positions:
                ref_pos = random.choice(list(existing_positions))
                pos = ref_pos + random.randint(-2, 2)
            else:
                pos = random.randint(50, len(wt_sequence) - 50)
            
            pos_idx = pos - 1
            
            if pos_idx < 0 or pos_idx >= len(variant_seq):
                continue
            
            if pos in self.functional_exclusion:
                continue
            
            wt_aa = wt_sequence[pos_idx]
            current_aa = variant_seq[pos_idx]
            
            if wt_aa != current_aa:
                continue
            
            if wt_aa in stabilizing:
                to_aa = random.choice(stabilizing[wt_aa])
                
                in_tm = any(pos in tm for tm in self.membrane_regions.values())
                if in_tm and to_aa == 'P':
                    continue
                
                variant_seq[pos_idx] = to_aa
                synergistic.append(f'{wt_aa}{pos}{to_aa}')
        
        return synergistic
    
    def _heuristic_fallback(self, wt_sequence: str, n_variants: int) -> List[Dict]:
        """Pure heuristic variants when no Track 1 data"""
        
        variants = []
        
        helix_stabilizers = {
            'G': 'A', 'S': 'T', 'V': 'I', 'L': 'I'
        }
        
        for i in range(n_variants):
            variant_seq = list(wt_sequence)
            variant_muts = []
            
            n_muts = np.random.randint(3, 6)
            
            for _ in range(n_muts):
                tm_positions = []
                for tm_range in self.membrane_regions.values():
                    tm_positions.extend(list(tm_range))
                
                if tm_positions:
                    pos = random.choice(tm_positions)
                    pos_idx = pos - 1
                    
                    if pos in self.functional_exclusion:
                        continue
                    
                    if pos_idx >= len(variant_seq):
                        continue
                    
                    wt_aa = variant_seq[pos_idx]
                    
                    if wt_aa in helix_stabilizers:
                        to_aa = helix_stabilizers[wt_aa]
                        variant_seq[pos_idx] = to_aa
                        variant_muts.append(f'{wt_aa}{pos}{to_aa}')
            
            variants.append({
                'id': f'T2_HEURISTIC_{i:04d}',
                'sequence': ''.join(variant_seq),
                'mutations': variant_muts,
                'method': 'heuristic',
                'based_on_track1': False
            })
        
        return variants
    
    async def _multi_objective_score(
        self, variant: Dict, wt_sequence: str, target_temp: int
    ) -> Dict:
        """Multi-objective scoring"""
        
        scores = {}
        
        # 1. Predict ΔTm
        if self.client:
            try:
                tm_result = await self.client.predict_tm(variant['sequence'], direction='Increase')
                predicted_tm = tm_result.get('tm', 42)
                scores['predicted_tm'] = predicted_tm
                scores['delta_tm'] = predicted_tm - 42
            except:
                scores['delta_tm'] = self._heuristic_delta_tm(variant)
                scores['predicted_tm'] = 42 + scores['delta_tm']
        else:
            scores['delta_tm'] = self._heuristic_delta_tm(variant)
            scores['predicted_tm'] = 42 + scores['delta_tm']
        
        # 2. ΔΔG_membrane
        scores['delta_dg_membrane'] = self._estimate_membrane_ddg(variant, wt_sequence)
        
        # 3. Helix propensity
        scores['helix_score'] = self._helix_propensity_score(variant, wt_sequence)
        
        # 4. Distance to functional sites
        scores['min_distance_functional'] = self._min_distance_to_functional(variant)
        
        # 5. Composite score
        scores['composite_score'] = (
            scores['delta_tm'] * 0.35 +
            -scores['delta_dg_membrane'] * 0.25 +
            scores['helix_score'] * 0.20 +
            (scores['min_distance_functional'] / 10) * 0.20
        )
        
        return scores
    
    def _heuristic_delta_tm(self, variant: Dict) -> float:
        """Estimate ΔTm from mutation patterns"""
        
        delta = 0.0
        mutations = variant.get('mutations', [])
        
        for mut in mutations:
            if len(mut) < 3:
                continue
            
            from_aa = mut[0]
            to_aa = mut[-1]
            
            if from_aa in 'GSN' and to_aa in 'AILV':
                delta += 1.5
            if from_aa == 'G' and to_aa == 'A':
                delta += 1.0
            if to_aa == 'P':
                delta += 0.5
        
        return min(delta, 8.0)
    
    def _estimate_membrane_ddg(self, variant: Dict, wt_sequence: str) -> float:
        """Approximate ΔΔG in membrane context"""
        
        ddg = 0.0
        mutations = variant.get('mutations', [])
        
        for mut in mutations:
            if len(mut) < 3:
                continue
            
            try:
                pos = int(mut[1:-1])
                in_tm = any(pos in tm for tm in self.membrane_regions.values())
                
                from_aa = mut[0]
                to_aa = mut[-1]
                
                if in_tm:
                    if from_aa in 'STNQ' and to_aa in 'AILMFWV':
                        ddg -= 1.0
                    elif from_aa in 'AILMFWV' and to_aa in 'STNQDE':
                        ddg += 1.5
            except:
                continue
        
        return ddg
    
    def _helix_propensity_score(self, variant: Dict, wt_sequence: str) -> float:
        """Score helix propensity changes"""
        
        helix_prop = {
            'A': 1.41, 'L': 1.34, 'M': 1.30, 'E': 1.26, 'K': 1.23,
            'F': 1.16, 'Q': 1.13, 'I': 1.09, 'W': 1.09, 'D': 1.04,
            'V': 0.98, 'Y': 0.74, 'R': 0.79, 'T': 0.76, 'S': 0.57,
            'C': 0.66, 'N': 0.76, 'H': 0.89, 'G': 0.43, 'P': 0.34
        }
        
        score_change = 0.0
        mutations = variant.get('mutations', [])
        
        for mut in mutations:
            if len(mut) < 3:
                continue
            
            try:
                pos = int(mut[1:-1])
                in_tm = any(pos in tm for tm in self.membrane_regions.values())
                
                if in_tm:
                    from_aa = mut[0]
                    to_aa = mut[-1]
                    
                    from_prop = helix_prop.get(from_aa, 1.0)
                    to_prop = helix_prop.get(to_aa, 1.0)
                    
                    score_change += (to_prop - from_prop)
            except:
                continue
        
        return score_change
    
    def _min_distance_to_functional(self, variant: Dict) -> int:
        """Minimum distance to any functional site"""
        
        mutations = variant.get('mutations', [])
        min_dist = 100
        
        for mut in mutations:
            if len(mut) < 3:
                continue
            
            try:
                pos = int(mut[1:-1])
                
                for func_pos in self.functional_exclusion:
                    dist = abs(pos - func_pos)
                    min_dist = min(min_dist, dist)
            except:
                continue
        
        return min_dist
    
    def _apply_mutations_to_sequence(self, sequence: str, mutations: List[str]) -> str:
        """Apply a list of mutations to a sequence"""
        
        seq = list(sequence)
        
        for mut in mutations:
            if len(mut) < 3:
                continue
            
            try:
                pos = int(mut[1:-1]) - 1
                to_aa = mut[-1]
                
                if 0 <= pos < len(seq):
                    seq[pos] = to_aa
            except:
                continue
        
        return ''.join(seq)
    
    def _get_mutations(self, wt: str, variant: str) -> List[str]:
        """Extract mutations between WT and variant"""
        
        mutations = []
        min_len = min(len(wt), len(variant))
        
        for i in range(min_len):
            if wt[i] != variant[i]:
                mutations.append(f'{wt[i]}{i+1}{variant[i]}')
        
        return mutations