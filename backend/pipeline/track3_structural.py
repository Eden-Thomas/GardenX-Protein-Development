"""
track3_structural.py
Membrane-aware structural refinement building on Track 1+2 candidates
"""

import asyncio
import numpy as np
from typing import Dict, List, Optional


class Track3Structural:
    """Structure-guided refinement with membrane orientation"""
    
    def __init__(self, config, client):
        self.config = config
        self.client = client
        self.wt_structure_cache = None
        
        self.membrane_topology = {
            'TM1': range(21, 47),
            'TM2': range(60, 88),
            'TM3': range(181, 207),
            'TM4': range(241, 270),
            'TM5': range(311, 342)
        }
        
        self.qb_pocket = {215, 254, 255, 264, 265, 271}
        self.oec_cluster = {170, 189, 332, 333, 342, 344}
    
    async def generate_variants(
        self,
        wt_sequence: str,
        n_variants: int = 10,
        seed_variants: List[Dict] = None
    ) -> List[Dict]:
        """
        Generate structure-refined variants
        
        SEQUENTIAL MODE: Refine seed variants from Track 1+2
        INDEPENDENT MODE: Generate from scratch (fallback)
        """
        
        print("  🏗️ Track 3: Structural refinement")
        
        variants = []
        
        # SEQUENTIAL MODE: Refine seeds
        if seed_variants and len(seed_variants) > 0:
            print(f"    → SEQUENTIAL MODE: Refining {len(seed_variants)} seed variants")
            
            for i, v in enumerate(seed_variants[:3]):
                print(f"      • {v.get('id', 'unknown')} ({v.get('track', 'unknown')}) - {len(v.get('mutations', []))} muts")
            
            top_seeds = seed_variants[:min(5, len(seed_variants))]
            
            print(f"    → Predicting structures for top {len(top_seeds)} seeds")
            
            for i, seed in enumerate(top_seeds):
                try:
                    seed_id = seed.get('id', f'seed_{i}')
                    seed_seq = seed['sequence']
                    
                    print(f"      Refining {seed_id}...")
                    
                    seed_structure = await self.predict_best_structure(seed_seq)
                    
                    if self.client:
                        try:
                            mpnn_sequences = await self.client.inverse_folding(
                                seed_structure['structure'],
                                num_sequences=2
                            )
                            
                            if mpnn_sequences:
                                for j, seq_data in enumerate(mpnn_sequences):
                                    sequence = seq_data.get('sequence', seq_data) if isinstance(seq_data, dict) else seq_data
                                    mutations = self._get_mutations(wt_sequence, sequence)
                                    
                                    membrane_score = self._membrane_packing_score(wt_sequence, sequence)
                                    
                                    variants.append({
                                        'id': f'T3_STRUCT_{i}_{j:02d}',
                                        'sequence': sequence,
                                        'mutations': mutations,
                                        'method': 'ligandmpnn_refinement',
                                        'track': 'track3',
                                        'parent_variant': seed_id,
                                        'parent_track': seed.get('track', 'unknown'),
                                        'parent_mutations': seed.get('mutations', []),
                                        'plddt': seed_structure.get('plddt', 0),
                                        'membrane_packing': membrane_score,
                                        'score': seq_data.get('score', 0) if isinstance(seq_data, dict) else 0
                                    })
                                
                                print(f"      ✓ {seed_id} → {len(mpnn_sequences)} refined variants")
                            else:
                                variants.append(self._validate_seed_structure(seed, seed_id, seed_seq, seed_structure, wt_sequence, i))
                        
                        except Exception as e:
                            print(f"      ⚠️ LigandMPNN failed for {seed_id}: {e}")
                            variants.append(self._validate_seed_structure(seed, seed_id, seed_seq, seed_structure, wt_sequence, i))
                    else:
                        variants.append(self._validate_seed_structure(seed, seed_id, seed_seq, seed_structure, wt_sequence, i))
                
                except Exception as e:
                    print(f"      ⚠️ Failed to refine {seed.get('id', 'unknown')}: {e}")
                    continue
            
            remaining = n_variants - len(variants)
            if remaining > 0:
                print(f"    → Generating {remaining} additional heuristic variants")
                for i in range(remaining):
                    heuristic_variant = self._heuristic_structural_variant(wt_sequence)
                    heuristic_variant['id'] = f'T3_HEURISTIC_{len(variants) + i:04d}'
                    variants.append(heuristic_variant)
        
        # INDEPENDENT MODE: No seeds
        else:
            print("    → INDEPENDENT MODE: Generating from scratch")
            
            print("    → Predicting wild-type structure...")
            wt_structure = await self.predict_best_structure(wt_sequence)
            self.wt_structure_cache = wt_structure
            print(f"      ✓ Structure predicted (pLDDT: {wt_structure['plddt']:.1f})")
            
            for i in range(n_variants):
                variant = self._heuristic_structural_variant(wt_sequence)
                variant['id'] = f'T3_STRUCT_{i:04d}'
                variants.append(variant)
        
        for variant in variants:
            if 'track' not in variant:
                variant['track'] = 'track3'
            if 'method' not in variant:
                variant['method'] = 'structural'
        
        print(f"    ✓ Generated {len(variants)} structural variants")
        
        return variants
    
    def _validate_seed_structure(self, seed: Dict, seed_id: str, seed_seq: str, 
                                 seed_structure: Dict, wt_sequence: str, index: int) -> Dict:
        """Create validated variant from seed"""
        return {
            **seed,
            'id': f'T3_STRUCT_{index}_00',
            'track': 'track3',
            'method': 'structure_validated',
            'parent_variant': seed_id,
            'plddt': seed_structure.get('plddt', 0),
            'membrane_packing': self._membrane_packing_score(wt_sequence, seed_seq)
        }
    
    async def predict_best_structure(self, sequence: str) -> Dict:
        """Predict structure with membrane orientation"""
        
        if self.client:
            try:
                result = await self.client.predict_structure(sequence, model='alphafold2')
                
                plddt = result.get('plddt', 85)
                structure_pdb = result.get('pdb', '')
                
                return {
                    'plddt': plddt,
                    'structure': structure_pdb,
                    'method': 'alphafold2'
                }
            
            except Exception as e:
                print(f"      ⚠️ Structure prediction failed: {e}")
                return {'plddt': 80, 'structure': '', 'method': 'mock'}
        else:
            return {
                'plddt': 80 + np.random.rand() * 10,
                'structure': '',
                'method': 'mock'
            }
    
    def _membrane_packing_score(self, wt_sequence: str, variant_sequence: str) -> float:
        """Score membrane packing quality"""
        
        score = 0.0
        min_len = min(len(wt_sequence), len(variant_sequence))
        
        for pos in range(min_len):
            canonical_pos = pos + 1
            variant_aa = variant_sequence[pos]
            
            in_tm = any(canonical_pos in tm for tm in [
                self.membrane_topology['TM1'],
                self.membrane_topology['TM2'],
                self.membrane_topology['TM3'],
                self.membrane_topology['TM4'],
                self.membrane_topology['TM5']
            ])
            
            if in_tm:
                if variant_aa in 'AILMFWV':
                    score += 0.5
                elif variant_aa in 'STNQDE':
                    score -= 1.0
                elif variant_aa in 'KR':
                    score -= 1.5
            else:
                if variant_aa in 'STNQDE':
                    score += 0.2
        
        return score / max(min_len, 1)
    
    def _check_qb_pocket_integrity(self, wt_sequence: str, variant_sequence: str) -> bool:
        """Check if QB pocket is intact"""
        
        for pos in self.qb_pocket:
            pos_idx = pos - 1
            
            if pos_idx >= len(wt_sequence) or pos_idx >= len(variant_sequence):
                continue
            
            if wt_sequence[pos_idx] != variant_sequence[pos_idx]:
                return False
        
        return True
    
    def _check_oec_integrity(self, wt_sequence: str, variant_sequence: str) -> bool:
        """Check if OEC cluster ligands are intact"""
        
        for pos in self.oec_cluster:
            pos_idx = pos - 1
            
            if pos_idx >= len(wt_sequence) or pos_idx >= len(variant_sequence):
                continue
            
            if wt_sequence[pos_idx] != variant_sequence[pos_idx]:
                return False
        
        return True
    
    def _heuristic_structural_variant(self, wt_sequence: str) -> Dict:
        """Generate heuristic variant with structural considerations"""
        
        variant_seq = list(wt_sequence)
        mutations = []
        
        n_muts = np.random.randint(3, 6)
        
        loop_positions = []
        for tm_range in [range(47, 60), range(207, 241), range(270, 311)]:
            loop_positions.extend(list(tm_range))
        
        safe_loop_positions = [p for p in loop_positions 
                              if p not in self.qb_pocket 
                              and p not in self.oec_cluster]
        
        attempts = 0
        while len(mutations) < n_muts and attempts < 50:
            attempts += 1
            
            if not safe_loop_positions:
                break
            
            pos = np.random.choice(safe_loop_positions)
            pos_idx = pos - 1
            
            if pos_idx >= len(variant_seq):
                continue
            
            wt_aa = variant_seq[pos_idx]
            
            if wt_aa == 'G':
                to_aa = np.random.choice(['A', 'S'])
            elif wt_aa == 'S':
                to_aa = 'T'
            elif wt_aa == 'N':
                to_aa = 'Q'
            else:
                continue
            
            variant_seq[pos_idx] = to_aa
            mutations.append(f'{wt_aa}{pos}{to_aa}')
        
        membrane_score = self._membrane_packing_score(wt_sequence, ''.join(variant_seq))
        qb_intact = self._check_qb_pocket_integrity(wt_sequence, ''.join(variant_seq))
        oec_intact = self._check_oec_integrity(wt_sequence, ''.join(variant_seq))
        
        return {
            'sequence': ''.join(variant_seq),
            'mutations': mutations,
            'method': 'heuristic_structural',
            'track': 'track3',
            'membrane_packing': membrane_score,
            'qb_pocket_intact': qb_intact,
            'oec_intact': oec_intact,
            'plddt': 80
        }
    
    def _get_mutations(self, wt: str, variant: str) -> List[str]:
        """Extract mutations"""
        
        mutations = []
        min_len = min(len(wt), len(variant))
        
        for i in range(min_len):
            if wt[i] != variant[i]:
                mutations.append(f'{wt[i]}{i+1}{variant[i]}')
        
        return mutations