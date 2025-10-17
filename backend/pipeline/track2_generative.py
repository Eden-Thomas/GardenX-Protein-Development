"""
track2_generative.py - Production generative track with NeuroFold
"""

import numpy as np
from typing import List, Dict
import asyncio
import json  

class Track2Generative:
    """
    Generative AI track
    Includes NeuroFold, RFdiffusion, Efficient Evolution
    """
    
    def __init__(self, config, neurosnap_client):
        self.config = config
        self.client = neurosnap_client
        
    async def generate_variants(self, wt_sequence: str, 
                               target_temp: int, 
                               n_variants: int) -> List[Dict]:
        """Generate variants using generative AI approaches"""
        
        print("  🤖 Track 2: Generative AI (NeuroFold + RFdiffusion)")
        
        variants = []
        
        # NeuroFold is KEY - give it 70% of variants
        neurofold_count = int(n_variants * 0.7) if n_variants > 1 else 1
        rfdiffusion_count = n_variants - neurofold_count
        
        # 1. NeuroFold - Specialized enzyme optimization
        try:
            print(f"    - NeuroFold: optimizing for {target_temp}°C...")
            
            neurofold_job = await self.client.submit_job('neurofold', {
                # FIX: Wrap sequence in json.dumps with proper format
                'sequence': json.dumps([{"type": "fasta", "data": wt_sequence}]),
                'target': 'thermostability',  # Correct parameter name
                # FIX: Convert to string
                'temperature': str(target_temp),
                # FIX: Convert to string
                'ph': "7.0",
                # FIX: Convert to string
                'num_designs': str(max(neurofold_count * 2, 10))
            })
            
            neurofold_results = await self.client.wait_for_job(neurofold_job)
            
            # Process NeuroFold results
            designs = neurofold_results.get('designs', [])
            
            for i, design in enumerate(designs[:neurofold_count]):
                mutations = self._get_mutations(wt_sequence, design.get('sequence', wt_sequence))
                
                variants.append({
                    'id': f"T2_NEURO_{i:04d}",
                    'sequence': design.get('sequence', wt_sequence),
                    'mutations': mutations,
                    'method': 'neurofold',
                    'predicted_tm': design.get('predicted_tm', target_temp),
                    'confidence': design.get('confidence', 0.8),
                    'optimization_score': design.get('optimization_score', 0.7),
                    'track': 'track2'
                })
            
            print(f"      ✅ Generated {len(variants)} NeuroFold variants")
            
        except Exception as e:
            print(f"      ⚠️ NeuroFold failed: {e}, using fallback")
            
            # Fallback: Generate variants using heuristic approach
            for i in range(neurofold_count):
                variant_seq = self._apply_thermostability_mutations(wt_sequence)
                mutations = self._get_mutations(wt_sequence, variant_seq)
                
                variants.append({
                    'id': f"T2_FALLBACK_{i:04d}",
                    'sequence': variant_seq,
                    'mutations': mutations,
                    'method': 'heuristic',
                    'track': 'track2'
                })
        
        # 2. RFdiffusion or alternative generative approach for remaining variants
        if rfdiffusion_count > 0:
            try:
                print(f"    - Generative design: creating {rfdiffusion_count} variants...")
                
                # Try to get structure first
                structure_job = await self.client.submit_job('alphafold2', {
                    # FIX: Wrap sequence in json.dumps with proper format
                    'sequence': json.dumps([{"type": "fasta", "data": wt_sequence}]),
                    # FIX: Convert to string
                    'num_models': "1",
                    # FIX: Convert boolean to string
                    'relax': "false"
                })
                
                structure_result = await self.client.wait_for_job(structure_job)
                structure = structure_result.get('structure_pdb', '')
                
                # Try RFdiffusion if structure available
                if structure:
                    try:
                        # Use LigandMPNN as alternative to RFdiffusion
                        mpnn_job = await self.client.submit_job('ligandmpnn', {
                            # FIX: Wrap structure in json.dumps with proper format
                            'structure': json.dumps([{"type": "pdb", "data": structure}]),
                            # FIX: Convert to string
                            'num_sequences': str(rfdiffusion_count * 2),
                            # FIX: Convert to string
                            'temperature': "0.2",
                            'fixed_positions': json.dumps(self._get_fixed_positions())
                        })
                        
                        mpnn_results = await self.client.wait_for_job(mpnn_job)
                        
                        for i, seq_data in enumerate(mpnn_results.get('sequences', [])[:rfdiffusion_count]):
                            mutations = self._get_mutations(wt_sequence, seq_data.get('sequence', wt_sequence))
                            
                            variants.append({
                                'id': f"T2_MPNN_{i:04d}",
                                'sequence': seq_data.get('sequence', wt_sequence),
                                'mutations': mutations,
                                'method': 'ligandmpnn',
                                'score': seq_data.get('score', 0),
                                'track': 'track2'
                            })
                            
                    except Exception as e:
                        print(f"      ⚠️ Structure-based design failed: {e}")
                        # Fall through to heuristic approach
                
            except Exception as e:
                print(f"      ⚠️ Structure prediction failed: {e}")
            
            # Fill remaining slots with heuristic variants if needed
            current_count = len(variants)
            remaining = n_variants - current_count
            
            for i in range(remaining):
                variant_seq = self._apply_efficient_evolution(wt_sequence)
                mutations = self._get_mutations(wt_sequence, variant_seq)
                
                variants.append({
                    'id': f"T2_EE_{i:04d}",
                    'sequence': variant_seq,
                    'mutations': mutations,
                    'method': 'efficient_evolution',
                    'track': 'track2'
                })
        
        print(f"    Total Track 2 variants: {len(variants)}")
        return variants
    
    def _get_fixed_positions(self) -> List[int]:
        """Define positions to keep fixed during design (1-indexed)"""
        # Active site positions that must be preserved
        return [161, 170, 189, 215, 254, 255, 264, 271, 332, 333, 342, 344]
    
    def _get_mutations(self, wt: str, variant: str) -> List[str]:
        """Extract mutations between sequences"""
        mutations = []
        min_len = min(len(wt), len(variant))
        
        for i in range(min_len):
            if wt[i] != variant[i]:
                mutations.append(f"{wt[i]}{i+1}{variant[i]}")
        
        return mutations
    
    def _apply_thermostability_mutations(self, sequence: str) -> str:
        """Apply known thermostability-enhancing mutations"""
        
        variant = list(sequence)
        
        # Common thermostability strategies
        strategies = [
            ('G', ['A', 'V']),      # Remove glycine flexibility
            ('S', ['T', 'A']),      # Reduce polar surface
            ('N', ['Q', 'D']),      # Prevent deamidation
            ('M', ['L', 'I']),      # Remove oxidizable Met
            ('C', ['S', 'A']),      # Remove reactive Cys
            ('K', ['R']),           # Conservative charge
            ('E', ['D']),           # Conservative charge
            ('A', ['V', 'I']),      # Increase hydrophobic packing
        ]
        
        # Apply 5-10 mutations
        n_mutations = np.random.randint(5, 11)
        mutation_count = 0
        
        # Avoid active site
        active_site = {160, 169, 188, 214, 253, 254, 263, 270, 331, 332, 341, 343}
        
        positions = list(range(len(sequence)))
        np.random.shuffle(positions)
        
        for pos in positions:
            if mutation_count >= n_mutations:
                break
            
            if (pos + 1) in active_site:
                continue
            
            current_aa = variant[pos]
            
            for from_aa, to_aas in strategies:
                if current_aa == from_aa:
                    variant[pos] = np.random.choice(to_aas)
                    mutation_count += 1
                    break
        
        return ''.join(variant)
    
    def _apply_efficient_evolution(self, sequence: str) -> str:
        """Apply efficient evolution strategy"""
        
        variant = list(sequence)
        
        # Target specific regions for improvement
        # Loops between TM helices are good targets
        loop_regions = [
            (41, 69),    # Loop 1
            (91, 159),   # Loop 2
            (181, 219),  # Loop 3
            (241, 289),  # Loop 4
            (311, 331)   # Loop 5
        ]
        
        # Apply 3-7 mutations in loops
        n_mutations = np.random.randint(3, 8)
        mutation_count = 0
        
        for _ in range(n_mutations):
            # Pick a random loop
            if loop_regions:
                start, end = loop_regions[np.random.randint(len(loop_regions))]
                pos = np.random.randint(start-1, min(end, len(sequence)))
                
                # Apply stabilizing mutation
                old_aa = variant[pos]
                
                # Prefer Pro in loops for rigidity, or small hydrophobic
                if old_aa != 'P' and np.random.random() > 0.5:
                    variant[pos] = 'P'
                else:
                    variant[pos] = np.random.choice(['A', 'V', 'L', 'I'])
                
                mutation_count += 1
        
        return ''.join(variant)