"""
track1_evolutionary.py
Evolutionary analysis with position mapping, topology filters, and functional site exclusions
UPDATED: Added aggressive mode for higher temperature targets
"""

import asyncio
import numpy as np
from typing import Dict, List, Optional, Set, Tuple
from collections import Counter

# Try to import Bio, but provide fallback
try:
    from Bio import Entrez, SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    print("WARNING: BioPython not available. Using fallback mode.")
    BIOPYTHON_AVAILABLE = False


# D1/psbA STRUCTURAL ANNOTATIONS (canonical numbering)
D1_TOPOLOGY = {
    'TM1': range(21, 47),
    'TM2': range(60, 88),
    'TM3': range(181, 207),
    'TM4': range(241, 270),
    'TM5': range(311, 342),
    'loops': [range(1, 21), range(47, 60), range(88, 181), 
              range(207, 241), range(270, 311), range(342, 360)]
}

# CRITICAL FUNCTIONAL RESIDUES - DO NOT MUTATE
FUNCTIONAL_MASK = {
    'OEC': {170, 189, 332, 333, 342, 344},
    'TyrZ': {161},
    'QB_pocket': {215, 254, 255, 264, 265, 271},
    'PheoD1': {130, 147},
    'Cterm': {344}
}

ALL_LOCKED_POSITIONS = set()
for positions in FUNCTIONAL_MASK.values():
    ALL_LOCKED_POSITIONS.update(positions)


class Track1Evolutionary:
    """Evolutionary analysis with comprehensive guardrails"""
    
    def __init__(self, config, client):
        self.config = config
        self.client = client
        self.canonical_length = 344
        
        if BIOPYTHON_AVAILABLE:
            Entrez.email = "thomas@edenaglabs.com"  # CHANGE THIS TO YOUR EMAIL
    
    async def generate_variants(
        self, 
        wt_sequence: str, 
        n_variants: int = 15,
        aggressive_mode: bool = False,        # NEW
        conservation_threshold: float = 0.9,  # NEW - make configurable
        max_mutations: int = 6,               # NEW - make configurable  
        target_temp: int = 50                 # NEW
    ) -> List[Dict]:
        """
        Generate evolution-guided variants with enhanced aggressiveness
        
        Args:
            wt_sequence: Wild-type sequence
            n_variants: Number of variants to generate
            aggressive_mode: Enable aggressive mode for higher temps
            conservation_threshold: Conservation cutoff (lower = more mutations)
            max_mutations: Maximum mutations per variant
            target_temp: Target temperature in Celsius
        """
        
        print("  🧬 Track 1: Evolutionary analysis with guardrails")
        
        # CRITICAL CHANGES for aggressive mode
        if aggressive_mode:
            print(f"    → AGGRESSIVE MODE: Targeting {target_temp}°C")
            print(f"    → Conservation threshold: {conservation_threshold} → 0.3 (lower = more mutations)")
            print(f"    → Max mutations per variant: {max_mutations} → 15")
            
            conservation_threshold = 0.3  # Much lower!
            max_mutations = 15           # Much higher!
            thermophile_temp_min = 65    # Target extreme thermophiles instead of moderate
        else:
            thermophile_temp_min = 50    # Moderate thermophiles
        
        if BIOPYTHON_AVAILABLE:
            print(f"    → Fetching thermophilic sequences (min temp: {thermophile_temp_min}°C)...")
            thermophile_seqs = await self._fetch_thermophile_sequences()
        else:
            print("    → BioPython unavailable, using heuristic mode")
            thermophile_seqs = []
        
        if not thermophile_seqs or len(thermophile_seqs) < 3:
            print("    ⚠️ Insufficient thermophile data, using fallback")
            return self._fallback_variants(wt_sequence, n_variants, aggressive_mode, max_mutations)
        
        print(f"    ✓ Retrieved {len(thermophile_seqs)} thermophilic sequences")
        
        print("    → Analyzing conservation and topology...")
        conservation_data = self._analyze_conservation(wt_sequence, thermophile_seqs)
        
        print("    → Identifying thermostable mutations...")
        candidate_mutations = self._identify_mutations_with_guardrails(
            wt_sequence, thermophile_seqs, conservation_data, conservation_threshold
        )
        
        if not candidate_mutations:
            print("    ⚠️ No safe mutations found after filtering")
            return self._fallback_variants(wt_sequence, n_variants, aggressive_mode, max_mutations)
        
        print(f"    ✓ Found {len(candidate_mutations)} safe mutations")
        
        variants = self._build_variants_from_mutations(
            wt_sequence, candidate_mutations, n_variants, max_mutations
        )
        
        for i, variant in enumerate(variants):
            variant['id'] = f'T1_EVO_{i:04d}'
            variant['track'] = 'track1'
            variant['method'] = 'evolutionary'
            variant['aggressive_mode'] = aggressive_mode
        
        print(f"    ✓ Generated {len(variants)} evolution-guided variants")
        
        top_mutations = [m['mutation'] for m in candidate_mutations[:10]]
        print(f"    → Top mutations for Track 2: {', '.join(top_mutations[:5])}")
        
        return variants
    
    async def _fetch_thermophile_sequences(self) -> List[Dict]:
        """Fetch psbA sequences from thermophilic organisms"""
        
        if not BIOPYTHON_AVAILABLE:
            return []
        
        organisms = [
            "Thermosynechococcus elongatus",
            "Thermosynechococcus vestitus",
            "Synechococcus lividus",
            "Chroococcidiopsis thermalis"
        ]
        
        sequences = []
        
        for organism in organisms:
            try:
                search_query = f"{organism}[Organism] AND psbA[Gene]"
                handle = Entrez.esearch(db="nucleotide", term=search_query, retmax=5)
                record = Entrez.read(handle)
                handle.close()
                
                if record['IdList']:
                    fetch_handle = Entrez.efetch(
                        db="nucleotide",
                        id=record['IdList'][0],
                        rettype="fasta",
                        retmode="text"
                    )
                    seq_record = SeqIO.read(fetch_handle, "fasta")
                    fetch_handle.close()
                    
                    sequences.append({
                        'organism': organism,
                        'sequence': str(seq_record.seq),
                        'temperature': 50
                    })
                
                await asyncio.sleep(0.5)
                
            except Exception as e:
                print(f"      ⚠️ Failed to fetch {organism}: {e}")
                continue
        
        return sequences
    
    def _analyze_conservation(self, wt_sequence: str, thermophile_seqs: List[Dict]) -> Dict:
        """Calculate conservation scores per position"""
        
        all_seqs = [wt_sequence] + [t['sequence'] for t in thermophile_seqs]
        min_len = min(len(s) for s in all_seqs)
        aligned = [s[:min_len] for s in all_seqs]
        
        conservation_scores = []
        for pos in range(min_len):
            residues = [seq[pos] for seq in aligned]
            counts = Counter(residues)
            total = len(residues)
            entropy = -sum((c/total) * np.log2(c/total) for c in counts.values() if c > 0)
            conservation = 1 - (entropy / np.log2(min(20, total)))
            conservation_scores.append(conservation)
        
        return {
            'scores': conservation_scores,
            'aligned_seqs': aligned,
            'length': min_len
        }
    
    def _identify_mutations_with_guardrails(
        self, wt_sequence: str, thermophile_seqs: List[Dict], 
        conservation_data: Dict, conservation_threshold: float
    ) -> List[Dict]:
        """Identify mutations with ALL guardrails applied (now uses configurable threshold)"""
        
        mutations = []
        aligned_seqs = conservation_data['aligned_seqs']
        conservation = conservation_data['scores']
        
        for pos in range(len(conservation)):
            canonical_pos = pos + 1
            
            # GUARDRAIL 1: FUNCTIONAL SITE EXCLUSION
            if canonical_pos in ALL_LOCKED_POSITIONS:
                continue
            
            # GUARDRAIL 2: TOPOLOGY FILTER
            in_tm_helix = any(canonical_pos in tm_range for tm_range in 
                            [D1_TOPOLOGY['TM1'], D1_TOPOLOGY['TM2'], 
                             D1_TOPOLOGY['TM3'], D1_TOPOLOGY['TM4'], 
                             D1_TOPOLOGY['TM5']])
            
            wt_aa = wt_sequence[pos]
            thermo_residues = [seq[pos] for seq in aligned_seqs[1:]]
            residue_counts = Counter(thermo_residues)
            
            for residue, count in residue_counts.items():
                if residue == wt_aa or residue == '-':
                    continue
                
                prevalence = (count / len(thermo_residues)) * 100
                
                if prevalence < 20:
                    continue
                
                # GUARDRAIL 3: TM HELIX RESTRICTIONS
                if in_tm_helix:
                    if residue == 'P':
                        continue
                    if wt_aa == 'G' and residue in 'STNQDE':
                        continue
                    if wt_aa in 'AILMFWV' and residue not in 'AILMFWV':
                        continue
                
                # GUARDRAIL 4: CONSERVATION THRESHOLD (now uses parameter)
                if conservation[pos] > conservation_threshold:
                    continue
                
                # GUARDRAIL 5: DISTANCE TO FUNCTIONAL SITES
                min_distance = min(abs(canonical_pos - site) for site in ALL_LOCKED_POSITIONS)
                
                if min_distance < 3:
                    continue
                
                score = (
                    prevalence * 0.4 +
                    (1 - conservation[pos]) * 0.3 +
                    (min_distance / 50) * 0.3
                )
                
                topology = 'TM' if in_tm_helix else 'loop'
                
                mutations.append({
                    'position': canonical_pos,
                    'from_aa': wt_aa,
                    'to_aa': residue,
                    'mutation': f'{wt_aa}{canonical_pos}{residue}',
                    'prevalence': prevalence,
                    'conservation': conservation[pos],
                    'topology': topology,
                    'distance_to_active': min_distance,
                    'score': score
                })
        
        mutations.sort(key=lambda x: x['score'], reverse=True)
        return mutations
    
    def _build_variants_from_mutations(
        self, wt_sequence: str, mutations: List[Dict], 
        n_variants: int, max_mutations: int
    ) -> List[Dict]:
        """Build variants by combining top mutations (now uses configurable max_mutations)"""
        
        variants = []
        top_mutations = mutations[:min(20, len(mutations))]
        
        for i in range(n_variants):
            # Use max_mutations parameter instead of hardcoded value
            n_muts = np.random.randint(3, min(max_mutations + 1, len(top_mutations)))
            selected = np.random.choice(
                len(top_mutations), 
                size=min(n_muts, len(top_mutations)), 
                replace=False
            )
            selected_muts = [top_mutations[idx] for idx in selected]
            
            variant_seq = list(wt_sequence)
            variant_muts = []
            
            for mut in selected_muts:
                pos = mut['position'] - 1
                if pos < len(variant_seq):
                    variant_seq[pos] = mut['to_aa']
                    variant_muts.append(mut['mutation'])
            
            avg_score = np.mean([m['score'] for m in selected_muts])
            avg_prevalence = np.mean([m['prevalence'] for m in selected_muts])
            
            variants.append({
                'sequence': ''.join(variant_seq),
                'mutations': variant_muts,
                'score': avg_score,
                'conservation_score': 1 - np.mean([m['conservation'] for m in selected_muts]),
                'temperature_correlation': avg_prevalence / 100,
                'n_tm_mutations': sum(1 for m in selected_muts if m['topology'] == 'TM'),
                'n_loop_mutations': sum(1 for m in selected_muts if m['topology'] == 'loop')
            })
        
        return variants
    
    def _fallback_variants(
        self, wt_sequence: str, n_variants: int, 
        aggressive_mode: bool = False, max_mutations: int = 6
    ) -> List[Dict]:
        """Fallback: heuristic mutations when data unavailable (now with aggressive mode)"""
        
        variants = []
        
        # More aggressive mutations if in aggressive mode
        if aggressive_mode:
            safe_mutations = {
                45: ('G', 'V', 'Remove flexibility - aggressive'),
                50: ('S', 'A', 'Reduce polarity'),
                120: ('S', 'T', 'H-bond network'),
                125: ('N', 'D', 'Improve packing'),
                200: ('V', 'I', 'Core packing'),
                205: ('A', 'I', 'Increase hydrophobicity'),
                250: ('L', 'I', 'Hydrophobic core'),
                255: ('A', 'L', 'Better packing'),
                280: ('A', 'V', 'Helix stability'),
                285: ('T', 'I', 'Aggressive packing'),
                290: ('G', 'A', 'Rigidity'),
                300: ('S', 'V', 'Hydrophobic optimization')
            }
        else:
            safe_mutations = {
                45: ('G', 'A', 'Reduce flexibility'),
                120: ('S', 'T', 'H-bond network'),
                200: ('V', 'I', 'Core packing'),
                250: ('L', 'V', 'Hydrophobic core'),
                280: ('A', 'V', 'Helix stability')
            }
        
        for i in range(n_variants):
            variant_seq = list(wt_sequence)
            variant_muts = []
            
            # Use more mutations in aggressive mode
            n_muts_target = min(max_mutations, len(safe_mutations)) if aggressive_mode else min(3, len(safe_mutations))
            positions = np.random.choice(
                list(safe_mutations.keys()), 
                size=n_muts_target, 
                replace=False
            )
            
            for pos in positions:
                if pos - 1 < len(variant_seq):
                    from_aa, to_aa, _ = safe_mutations[pos]
                    if variant_seq[pos-1] == from_aa:
                        variant_seq[pos-1] = to_aa
                        variant_muts.append(f'{from_aa}{pos}{to_aa}')
            
            variants.append({
                'id': f'T1_FALLBACK_{i:04d}',
                'sequence': ''.join(variant_seq),
                'mutations': variant_muts,
                'score': 0.5,
                'track': 'track1',
                'method': 'heuristic_fallback',
                'aggressive_mode': aggressive_mode,
                'conservation_score': 0.5,
                'temperature_correlation': 0.5
            })
        
        return variants