"""
functional_design.py
Functional validation and combinatorial mutation design
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from itertools import combinations, product
from dataclasses import dataclass

@dataclass
class FunctionalScore:
    """Functional assessment of a variant"""
    photosynthesis_activity: float  # 0-1
    qb_binding_affinity: float  # 0-1
    electron_transfer_rate: float  # 0-1
    structural_integrity: float  # 0-1
    overall_function: float  # Combined score


class FunctionalValidator:
    """Validate that mutations preserve D1 function"""
    
    # Distance thresholds (Angstroms)
    QB_SITE_RADIUS = 12.0
    MN_CLUSTER_RADIUS = 10.0
    TYRZ_RADIUS = 8.0
    
    # Critical positions for D1 function
    CRITICAL_POSITIONS = {
        130, 147, 161, 170, 189, 215, 254, 255, 264, 271, 332, 333, 342, 344
    }
    
    def __init__(self):
        try:
            from structure_aware_design import CriticalResidueDatabase
            self.critical_db = CriticalResidueDatabase()
        except:
            # Fallback if structure module not available
            self.critical_db = None
    
    def predict_functional_impact(self,
                                  mutation: Dict,
                                  structural_context: Optional[Dict] = None) -> FunctionalScore:
        """
        Predict functional impact of mutation
        
        Args:
            mutation: Mutation dictionary with position, from, to
            structural_context: Optional structural information
        
        Returns:
            FunctionalScore object with detailed assessments
        """
        
        position = mutation.get('position', 0)
        from_aa = mutation.get('from', mutation.get('reference_aa', ''))
        to_aa = mutation.get('to', mutation.get('mutant_aa', ''))
        
        # Provide default structural context if missing
        if structural_context is None or not isinstance(structural_context, dict):
            structural_context = {}
        
        # Initialize scores
        qb_score = 1.0
        et_score = 1.0
        structure_score = 1.0
        
        # Check QB binding site impact
        # Use .get() with large default distance (far from active site)
        distance_to_active = structural_context.get('distance_to_active_site', 999.0)
        
        if distance_to_active < self.QB_SITE_RADIUS:
            # Near QB site - check if charge/polarity preserved
            charge_preserved = self._check_charge_preservation(from_aa, to_aa)
            if not charge_preserved:
                qb_score *= 0.6  # Significant impact
            else:
                qb_score *= 0.9  # Minor impact
        
        # Check electron transfer pathway
        if position in range(160, 170):  # Near TyrZ
            # Critical for electron transfer
            size_change = abs(self._get_residue_size(to_aa) - self._get_residue_size(from_aa))
            if size_change > 2:
                et_score *= 0.7
        
        # Check structural integrity
        in_membrane = structural_context.get('in_membrane', True)  # Assume membrane if unknown
        
        if in_membrane:
            # Membrane region - check hydrophobicity preservation
            hydrophobic_preserved = self._check_hydrophobicity_match(from_aa, to_aa)
            if not hydrophobic_preserved:
                structure_score *= 0.75
        
        # Check critical residues
        is_critical = False
        if self.critical_db and self.critical_db.is_critical(position):
            is_critical = True
        elif position in self.CRITICAL_POSITIONS:
            is_critical = True
        
        if is_critical:
            # Any mutation to critical residue is risky
            return FunctionalScore(
                photosynthesis_activity=0.0,
                qb_binding_affinity=0.0,
                electron_transfer_rate=0.0,
                structural_integrity=0.0,
                overall_function=0.0
            )
        
        # Calculate overall photosynthesis activity
        ps_activity = (qb_score + et_score + structure_score) / 3
        
        return FunctionalScore(
            photosynthesis_activity=ps_activity,
            qb_binding_affinity=qb_score,
            electron_transfer_rate=et_score,
            structural_integrity=structure_score,
            overall_function=ps_activity * structure_score
        )
    
    def _check_charge_preservation(self, aa1: str, aa2: str) -> bool:
        """Check if charge is preserved"""
        positive = set('KR')
        negative = set('DE')
        neutral = set('ACFGHILMNPQSTVWY')
        
        aa1_charge = 'pos' if aa1 in positive else 'neg' if aa1 in negative else 'neutral'
        aa2_charge = 'pos' if aa2 in positive else 'neg' if aa2 in negative else 'neutral'
        
        return aa1_charge == aa2_charge
    
    def _check_hydrophobicity_match(self, aa1: str, aa2: str) -> bool:
        """Check if hydrophobicity is similar"""
        hydrophobic = set('AILMFWV')
        return (aa1 in hydrophobic) == (aa2 in hydrophobic)
    
    def _get_residue_size(self, aa: str) -> int:
        """Get residue size scale"""
        size_scale = {
            'G': 1, 'A': 2, 'S': 2, 'C': 2, 'D': 3, 'P': 3, 'N': 3, 'T': 3,
            'E': 4, 'V': 4, 'Q': 4, 'H': 4, 'M': 4, 'I': 4, 'L': 4, 'K': 4,
            'R': 5, 'F': 5, 'Y': 5, 'W': 6
        }
        return size_scale.get(aa, 3)
    
    def filter_by_function(self,
                          mutations: List[Dict],
                          min_function_score: float = 0.7) -> List[Dict]:
        """
        Filter mutations to keep only functional ones
        
        Args:
            mutations: List of mutation dictionaries
            min_function_score: Minimum functional score to keep (0-1)
        
        Returns:
            List of functional mutations with scores added
        """
        
        functional_mutations = []
        
        for mutation in mutations:
            # Get structural context if available
            structural_context = mutation.get('structural_features', None)
            
            # If no structural context, provide safe defaults
            if structural_context is None or not isinstance(structural_context, dict):
                structural_context = {
                    'distance_to_active_site': 999.0,  # Far from active site
                    'in_membrane': True,  # Assume membrane region
                    'secondary_structure': 'unknown'
                }
            
            # Predict functional impact
            func_score = self.predict_functional_impact(mutation, structural_context)
            
            # Filter by minimum score
            if func_score.overall_function >= min_function_score:
                mutation['functional_score'] = func_score.overall_function
                mutation['function_details'] = {
                    'photosynthesis': func_score.photosynthesis_activity,
                    'qb_binding': func_score.qb_binding_affinity,
                    'electron_transfer': func_score.electron_transfer_rate,
                    'structural': func_score.structural_integrity
                }
                functional_mutations.append(mutation)
        
        return functional_mutations


class CombinatorialDesigner:
    """Design combinatorial mutation variants"""
    
    def __init__(self):
        self.validator = FunctionalValidator()
    
    def generate_combinations(self,
                             single_mutations: List[Dict],
                             max_mutations: int = 3,
                             min_spacing: int = 15) -> List[Dict]:
        """Generate combinations of mutations"""
        
        print(f"Generating combinations from {len(single_mutations)} mutations...")
        
        # Filter to top candidates
        top_mutations = sorted(
            single_mutations,
            key=lambda x: x.get('final_score', x.get('confidence', 0)),
            reverse=True
        )[:20]  # Top 20 only
        
        all_variants = []
        
        # Single mutations (baseline)
        for mut in top_mutations:
            mutation_id = mut.get('mutation', f"{mut.get('from', '?')}{mut.get('position', 0)}{mut.get('to', '?')}")
            all_variants.append({
                'mutations': [mut],
                'n_mutations': 1,
                'variant_id': f"single_{mutation_id}"
            })
        
        # Double combinations
        for mut1, mut2 in combinations(top_mutations, 2):
            # Check spatial spacing
            pos1 = mut1.get('position', 0)
            pos2 = mut2.get('position', 0)
            
            if abs(pos1 - pos2) < min_spacing:
                continue  # Too close - might interfere
            
            mut1_id = mut1.get('mutation', f"{mut1.get('from', '?')}{pos1}{mut1.get('to', '?')}")
            mut2_id = mut2.get('mutation', f"{mut2.get('from', '?')}{pos2}{mut2.get('to', '?')}")
            
            all_variants.append({
                'mutations': [mut1, mut2],
                'n_mutations': 2,
                'variant_id': f"double_{mut1_id}_{mut2_id}"
            })
        
        # Triple combinations (if requested)
        if max_mutations >= 3:
            for mut1, mut2, mut3 in combinations(top_mutations[:10], 3):  # Top 10 only for triples
                # Check pairwise spacing
                positions = [
                    mut1.get('position', 0),
                    mut2.get('position', 0),
                    mut3.get('position', 0)
                ]
                
                if not self._check_spacing(positions, min_spacing):
                    continue
                
                mut1_id = mut1.get('mutation', f"{mut1.get('from', '?')}{positions[0]}{mut1.get('to', '?')}")
                mut2_id = mut2.get('mutation', f"{mut2.get('from', '?')}{positions[1]}{mut2.get('to', '?')}")
                mut3_id = mut3.get('mutation', f"{mut3.get('from', '?')}{positions[2]}{mut3.get('to', '?')}")
                
                all_variants.append({
                    'mutations': [mut1, mut2, mut3],
                    'n_mutations': 3,
                    'variant_id': f"triple_{mut1_id}_{mut2_id}_{mut3_id}"
                })
        
        print(f"Generated {len(all_variants)} combinatorial variants")
        return all_variants
    
    def _check_spacing(self, positions: List[int], min_spacing: int) -> bool:
        """Check if all positions are sufficiently spaced"""
        for i in range(len(positions)):
            for j in range(i+1, len(positions)):
                if abs(positions[i] - positions[j]) < min_spacing:
                    return False
        return True
    
    def predict_combination_effect(self,
                                   variant: Dict,
                                   reference_sequence: str) -> Dict:
        """Predict effect of mutation combination"""
        
        mutations = variant['mutations']
        
        # Build mutant sequence
        mutant_seq = list(reference_sequence)
        for mut in mutations:
            pos_idx = mut.get('position', 0) - 1
            to_aa = mut.get('to', mut.get('mutant_aa', ''))
            
            if 0 <= pos_idx < len(mutant_seq) and to_aa:
                mutant_seq[pos_idx] = to_aa
        
        mutant_seq = ''.join(mutant_seq)
        
        # Predict cumulative effects
        cumulative_ddg = sum(m.get('delta_stability', 0) for m in mutations)
        
        # Check for synergy/antagonism
        # Simple model: assume slight positive synergy for well-spaced mutations
        if len(mutations) > 1:
            spacing_bonus = 0
            positions = [m.get('position', 0) for m in mutations]
            
            # Calculate average spacing
            spacings = [abs(positions[i] - positions[j]) 
                       for i in range(len(positions)) 
                       for j in range(i+1, len(positions))]
            
            if spacings:  # Avoid division by zero
                avg_spacing = np.mean(spacings)
                
                if avg_spacing > 30:  # Well-spaced
                    spacing_bonus = -0.3  # Additional stabilization
                
                cumulative_ddg += spacing_bonus
        
        # Predict functional retention
        min_function = min(m.get('functional_score', 0.9) for m in mutations)
        avg_function = np.mean([m.get('functional_score', 0.9) for m in mutations])
        
        # Combined function score (use minimum as bottleneck)
        combined_function = min_function * 0.7 + avg_function * 0.3
        
        return {
            'variant_sequence': mutant_seq,
            'predicted_ddg': cumulative_ddg,
            'estimated_tm_increase': -cumulative_ddg * 2.0,
            'functional_score': combined_function,
            'confidence': self._calculate_combination_confidence(variant),
            'recommendation': self._get_recommendation(cumulative_ddg, combined_function)
        }
    
    def _calculate_combination_confidence(self, variant: Dict) -> float:
        """Calculate confidence in combination prediction"""
        mutations = variant['mutations']
        
        # Average individual confidence
        individual_conf = np.mean([m.get('confidence', 0.5) for m in mutations])
        
        # Penalty for multiple mutations (less confident)
        n_penalty = 0.95 ** len(mutations)
        
        return individual_conf * n_penalty
    
    def _get_recommendation(self, ddg: float, function: float) -> str:
        """Get recommendation for variant"""
        if function < 0.6:
            return "HIGH_RISK_FUNCTION"
        elif ddg > 0:
            return "DESTABILIZING"
        elif ddg < -3.0 and function > 0.8:
            return "EXCELLENT_CANDIDATE"
        elif ddg < -1.5 and function > 0.7:
            return "GOOD_CANDIDATE"
        else:
            return "MODERATE_CANDIDATE"
    
    def rank_variants(self,
                     variants: List[Dict],
                     reference_sequence: str) -> List[Dict]:
        """Rank combinatorial variants"""
        
        scored_variants = []
        
        for variant in variants:
            prediction = self.predict_combination_effect(variant, reference_sequence)
            
            # Combined score: stability + function + confidence
            score = (
                max(0, -prediction['predicted_ddg'] / 5) * 0.4 +  # Stability (40%)
                prediction['functional_score'] * 0.4 +             # Function (40%)
                prediction['confidence'] * 0.2                     # Confidence (20%)
            )
            
            variant.update(prediction)
            variant['overall_score'] = score
            scored_variants.append(variant)
        
        # Sort by overall score
        scored_variants.sort(key=lambda x: x['overall_score'], reverse=True)
        
        return scored_variants


class ExperimentalPrioritizer:
    """Prioritize variants for experimental testing"""
    
    def __init__(self):
        self.designer = CombinatorialDesigner()
    
    def create_experimental_plan(self,
                                ranked_variants: List[Dict],
                                budget: int = 20) -> Dict:
        """Create experimental testing plan"""
        
        # Strategy: Test diverse set covering different approaches
        
        plan = {
            'wave_1_singles': [],      # 5 top single mutations
            'wave_2_doubles': [],      # 10 top double combinations
            'wave_3_triples': [],      # 5 top triple combinations
            'controls': []
        }
        
        # Wave 1: Top singles (build understanding)
        singles = [v for v in ranked_variants if v['n_mutations'] == 1]
        plan['wave_1_singles'] = singles[:5]
        
        # Wave 2: Top doubles (test synergy)
        doubles = [v for v in ranked_variants if v['n_mutations'] == 2]
        plan['wave_2_doubles'] = doubles[:10]
        
        # Wave 3: Top triples (push limits)
        triples = [v for v in ranked_variants if v['n_mutations'] == 3]
        plan['wave_3_triples'] = triples[:5]
        
        return plan
    
    def export_sequences_for_synthesis(self,
                                      experimental_plan: Dict,
                                      output_file: str = "synthesis_order.fasta"):
        """Export sequences in FASTA format for gene synthesis"""
        
        with open(output_file, 'w') as f:
            variant_num = 1
            
            for wave, variants in experimental_plan.items():
                if wave == 'controls':
                    continue
                
                for variant in variants:
                    variant_seq = variant.get('variant_sequence', '')
                    variant_id = variant.get('variant_id', f'variant{variant_num}')
                    
                    f.write(f">{wave}_variant{variant_num}_{variant_id}\n")
                    f.write(f"{variant_seq}\n")
                    variant_num += 1
        
        print(f"Exported {variant_num-1} sequences to {output_file}")
        return output_file


# Integration with main pipeline
def add_functional_and_combinatorial_design(
    filtered_mutations: List[Dict],
    reference_sequence: str,
    neurosnap_api
) -> Dict:
    """Complete functional validation and combinatorial design"""
    
    print("\n" + "="*60)
    print("FUNCTIONAL VALIDATION & COMBINATORIAL DESIGN")
    print("="*60)
    
    # Step 1: Functional filtering
    print("\n[1] Filtering by functional preservation...")
    validator = FunctionalValidator()
    functional_mutations = validator.filter_by_function(filtered_mutations)
    print(f"   Kept {len(functional_mutations)}/{len(filtered_mutations)} functional mutations")
    
    # Step 2: Generate combinations
    print("\n[2] Generating combinatorial variants...")
    designer = CombinatorialDesigner()
    variants = designer.generate_combinations(functional_mutations, max_mutations=3)
    
    # Step 3: Predict combination effects
    print("\n[3] Predicting combination effects...")
    ranked_variants = designer.rank_variants(variants, reference_sequence)
    
    # Step 4: Create experimental plan
    print("\n[4] Creating experimental testing plan...")
    prioritizer = ExperimentalPrioritizer()
    exp_plan = prioritizer.create_experimental_plan(ranked_variants, budget=20)
    
    # Step 5: Export for synthesis
    print("\n[5] Exporting sequences for gene synthesis...")
    fasta_file = prioritizer.export_sequences_for_synthesis(exp_plan)
    
    return {
        'functional_mutations': functional_mutations,
        'all_variants': ranked_variants,
        'experimental_plan': exp_plan,
        'synthesis_file': fasta_file,
        'summary': {
            'n_functional_singles': len([v for v in ranked_variants if v['n_mutations'] == 1]),
            'n_doubles': len([v for v in ranked_variants if v['n_mutations'] == 2]),
            'n_triples': len([v for v in ranked_variants if v['n_mutations'] == 3]),
            'top_variant': ranked_variants[0] if ranked_variants else None
        }
    }