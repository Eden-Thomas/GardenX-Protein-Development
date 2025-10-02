"""
structure_aware_design.py
Structure-aware mutation design using AlphaFold2 and spatial analysis
"""

import numpy as np
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import json

@dataclass
class StructuralContext:
    """Structural context for a residue"""
    position: int
    residue: str
    secondary_structure: str  # H=helix, E=sheet, C=coil
    solvent_accessibility: float  # 0-1
    distance_to_active_site: float  # Angstroms
    contacts: List[int]  # Nearby residues
    buried_surface_area: float  # Å²
    in_membrane: bool


class CriticalResidueDatabase:
    """Database of critical functional residues in D1"""
    
    # Based on literature and structural data
    CRITICAL_SITES = {
        'QB_binding': {
            'residues': [215, 252, 264, 266],  # His215, Phe252, Ser264, Phe266
            'description': 'Quinone QB binding pocket',
            'tolerance': 'zero'  # Never mutate
        },
        'metal_cluster': {
            'residues': [170, 189, 333, 344],  # Mn4CaO5 cluster ligands
            'description': 'Manganese cluster coordination',
            'tolerance': 'zero'
        },
        'tyrZ': {
            'residues': [161],  # Tyr161
            'description': 'TyrZ electron donor',
            'tolerance': 'zero'
        },
        'chlorophyll': {
            'residues': [197, 200],  # Chlorophyll binding
            'description': 'Chlorophyll coordination',
            'tolerance': 'low'  # Careful mutations only
        },
        'pheophytin': {
            'residues': [130, 133],
            'description': 'Pheophytin coordination',
            'tolerance': 'low'
        },
        'membrane_spanning': {
            'residues': list(range(20, 45)) + list(range(80, 105)) + 
                       list(range(150, 175)) + list(range(200, 225)) + 
                       list(range(270, 295)),
            'description': 'Transmembrane helices',
            'tolerance': 'medium'  # Hydrophobic substitutions OK
        }
    }
    
    @classmethod
    def is_critical(cls, position: int, site_type: str = 'any') -> bool:
        """Check if position is critical"""
        if site_type == 'any':
            for site_data in cls.CRITICAL_SITES.values():
                if position in site_data['residues']:
                    return True
            return False
        else:
            return position in cls.CRITICAL_SITES.get(site_type, {}).get('residues', [])
    
    @classmethod
    def get_tolerance(cls, position: int) -> str:
        """Get mutation tolerance for position"""
        for site_data in cls.CRITICAL_SITES.values():
            if position in site_data['residues']:
                return site_data['tolerance']
        return 'high'  # Not critical


class StructureAnalyzer:
    """Analyze protein structure from AlphaFold2 predictions"""
    
    def __init__(self):
        self.structure_cache = {}
    
    def parse_pdb_structure(self, pdb_content: str) -> Dict:
        """Parse PDB structure and extract key features"""
        
        # Parse ATOM records
        atoms = []
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                atom = {
                    'atom_num': int(line[6:11].strip()),
                    'atom_name': line[12:16].strip(),
                    'residue': line[17:20].strip(),
                    'chain': line[21],
                    'residue_num': int(line[22:26].strip()),
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'occupancy': float(line[54:60].strip()),
                    'b_factor': float(line[60:66].strip())  # pLDDT in AlphaFold
                }
                atoms.append(atom)
        
        return {'atoms': atoms}
    
    def calculate_distance_matrix(self, atoms: List[Dict]) -> np.ndarray:
        """Calculate Ca-Ca distance matrix"""
        # Get Ca atoms only
        ca_atoms = [a for a in atoms if a['atom_name'] == 'CA']
        
        n = len(ca_atoms)
        dist_matrix = np.zeros((n, n))
        
        for i in range(n):
            for j in range(i+1, n):
                pos_i = np.array([ca_atoms[i]['x'], ca_atoms[i]['y'], ca_atoms[i]['z']])
                pos_j = np.array([ca_atoms[j]['x'], ca_atoms[j]['y'], ca_atoms[j]['z']])
                
                dist = np.linalg.norm(pos_i - pos_j)
                dist_matrix[i, j] = dist
                dist_matrix[j, i] = dist
        
        return dist_matrix, ca_atoms
    
    def find_contacts(self, distance_matrix: np.ndarray, threshold: float = 8.0) -> Dict:
        """Find residue contacts within threshold distance"""
        contacts = {}
        
        for i in range(len(distance_matrix)):
            nearby = np.where((distance_matrix[i] < threshold) & (distance_matrix[i] > 0))[0]
            contacts[i] = nearby.tolist()
        
        return contacts
    
    def get_structural_context(self, position: int, structure_data: Dict) -> StructuralContext:
        """Get full structural context for a residue"""
        
        atoms = structure_data['atoms']
        dist_matrix, ca_atoms = self.calculate_distance_matrix(atoms)
        contacts = self.find_contacts(dist_matrix)
        
        # Find the Ca atom for this position
        residue_idx = position - 1  # Convert to 0-indexed
        
        if residue_idx >= len(ca_atoms):
            raise ValueError(f"Position {position} not found in structure")
        
        ca_atom = ca_atoms[residue_idx]
        
        # Calculate distance to critical sites
        qb_position = 215  # His215
        qb_idx = qb_position - 1
        dist_to_qb = dist_matrix[residue_idx, qb_idx] if qb_idx < len(ca_atoms) else 999
        
        # Predict secondary structure from B-factors and geometry
        # (Simplified - in real implementation use DSSP)
        ss = self._predict_secondary_structure(residue_idx, ca_atoms, dist_matrix)
        
        # Solvent accessibility (approximation from contacts)
        n_contacts = len(contacts.get(residue_idx, []))
        solvent_accessibility = max(0, 1 - n_contacts / 20)  # Rough estimate
        
        # Membrane detection (based on position in sequence)
        in_membrane = CriticalResidueDatabase.get_tolerance(position) == 'medium'
        
        return StructuralContext(
            position=position,
            residue=ca_atom['residue'],
            secondary_structure=ss,
            solvent_accessibility=solvent_accessibility,
            distance_to_active_site=dist_to_qb,
            contacts=contacts.get(residue_idx, []),
            buried_surface_area=n_contacts * 20,  # Rough estimate
            in_membrane=in_membrane
        )
    
    def _predict_secondary_structure(self, idx: int, ca_atoms: List, dist_matrix: np.ndarray) -> str:
        """Simple secondary structure prediction"""
        # Look at local geometry
        if idx < 2 or idx >= len(ca_atoms) - 2:
            return 'C'  # Coil at ends
        
        # Check if in regular pattern (simplified alpha helix detection)
        i_to_i4 = dist_matrix[idx, min(idx+4, len(ca_atoms)-1)]
        
        if 5.0 < i_to_i4 < 6.5:  # Typical alpha helix
            return 'H'
        elif i_to_i4 > 12:  # Extended
            return 'E'
        else:
            return 'C'


class StructureGuidedMutationDesigner:
    """Design mutations using structural context"""
    
    def __init__(self):
        self.analyzer = StructureAnalyzer()
        self.critical_db = CriticalResidueDatabase()
    
    def evaluate_mutation_safety(self, 
                                 position: int,
                                 from_aa: str,
                                 to_aa: str,
                                 structure_data: Optional[Dict] = None) -> Dict:
        """Evaluate if mutation is structurally safe"""
        
        # Rule 1: Critical residue check
        if self.critical_db.is_critical(position):
            tolerance = self.critical_db.get_tolerance(position)
            
            if tolerance == 'zero':
                return {
                    'safe': False,
                    'reason': 'Critical functional residue - NEVER mutate',
                    'risk_level': 'CRITICAL',
                    'recommendation': 'REJECT'
                }
        
        # Rule 2: Structure-based evaluation
        if structure_data:
            context = self.analyzer.get_structural_context(position, structure_data)
            
            # Check distance to active site
            if context.distance_to_active_site < 10.0:
                return {
                    'safe': False,
                    'reason': f'Too close to QB binding site ({context.distance_to_active_site:.1f} Å)',
                    'risk_level': 'HIGH',
                    'recommendation': 'CAUTION',
                    'context': context
                }
            
            # Buried residue check
            if context.solvent_accessibility < 0.2:
                # Buried - need similar properties
                if not self._check_property_match(from_aa, to_aa, 'hydrophobicity'):
                    return {
                        'safe': False,
                        'reason': 'Buried residue - hydrophobicity mismatch',
                        'risk_level': 'MEDIUM',
                        'recommendation': 'MODIFY',
                        'context': context
                    }
            
            # Secondary structure constraints
            if context.secondary_structure == 'H':  # Helix
                if to_aa == 'P':  # Proline breaks helices
                    return {
                        'safe': False,
                        'reason': 'Proline would disrupt helix',
                        'risk_level': 'MEDIUM',
                        'recommendation': 'REJECT',
                        'context': context
                    }
        
        # Passed all checks
        return {
            'safe': True,
            'reason': 'No structural conflicts detected',
            'risk_level': 'LOW',
            'recommendation': 'APPROVE',
            'context': context if structure_data else None
        }
    
    def _check_property_match(self, aa1: str, aa2: str, property_type: str) -> bool:
        """Check if amino acids match on a property"""
        
        properties = {
            'hydrophobicity': {
                'hydrophobic': set('AILMFWV'),
                'polar': set('STNQ'),
                'charged': set('KRED')
            }
        }
        
        if property_type not in properties:
            return True
        
        for group_name, group_set in properties[property_type].items():
            if (aa1 in group_set) == (aa2 in group_set):
                return True
        
        return False
    
    def filter_mutations_by_structure(self,
                                     mutations: List[Dict],
                                     structure_data: Dict) -> List[Dict]:
        """Filter mutation list based on structural safety"""
        
        filtered = []
        
        for mutation in mutations:
            safety = self.evaluate_mutation_safety(
                mutation['position'],
                mutation['from'],
                mutation['to'],
                structure_data
            )
            
            mutation['structure_safety'] = safety
            
            # Only keep safe or caution-level mutations
            if safety['recommendation'] in ['APPROVE', 'CAUTION']:
                mutation['final_score'] = mutation.get('score', 0.5) * (
                    1.0 if safety['recommendation'] == 'APPROVE' else 0.7
                )
                filtered.append(mutation)
        
        # Re-sort by final score
        filtered.sort(key=lambda x: x['final_score'], reverse=True)
        
        return filtered
    
    def generate_structure_aware_recommendations(self,
                                                sequence: str,
                                                neurosnap_api) -> List[Dict]:
        """Generate mutation recommendations using structure"""
        
        # Step 1: Get AlphaFold2 structure
        print("Predicting structure with AlphaFold2...")
        structure_result = neurosnap_api.predict_structure(sequence)
        pdb_content = structure_result['pdb_structure']
        structure_data = self.analyzer.parse_pdb_structure(pdb_content)
        
        # Step 2: Generate initial candidates (from comparative analysis)
        from analysis_engine import ComparativeAnalyzer
        analyzer = ComparativeAnalyzer()
        
        # This would connect to your existing analysis
        initial_candidates = []  # From your comparative analysis
        
        # Step 3: Filter by structural safety
        filtered_mutations = self.filter_mutations_by_structure(
            initial_candidates,
            structure_data
        )
        
        # Step 4: Add structural features to each
        for mutation in filtered_mutations:
            context = self.analyzer.get_structural_context(
                mutation['position'],
                structure_data
            )
            
            mutation['structural_features'] = {
                'secondary_structure': context.secondary_structure,
                'solvent_accessibility': context.solvent_accessibility,
                'distance_to_active_site': context.distance_to_active_site,
                'n_contacts': len(context.contacts),
                'in_membrane': context.in_membrane
            }
        
        return filtered_mutations


# Integrate with existing pipeline
def add_structure_awareness_to_pipeline(recommendations: List[Dict],
                                       sequence: str,
                                       neurosnap_api) -> List[Dict]:
    """Add structure-aware filtering to existing recommendations"""
    
    designer = StructureGuidedMutationDesigner()
    
    # Get structure
    print("Analyzing structure...")
    structure_result = neurosnap_api.predict_structure(sequence)
    pdb_content = structure_result['pdb_structure']
    structure_data = designer.analyzer.parse_pdb_structure(pdb_content)
    
    # Filter each recommendation
    filtered_recommendations = designer.filter_mutations_by_structure(
        recommendations,
        structure_data
    )
    
    print(f"Filtered {len(recommendations)} → {len(filtered_recommendations)} mutations")
    print(f"Rejected {len(recommendations) - len(filtered_recommendations)} due to structural conflicts")
    
    return filtered_recommendations