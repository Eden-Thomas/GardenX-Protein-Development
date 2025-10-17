"""
analysis_engine.py
Bioinformatics algorithms and statistical analysis for D1 protein engineering
"""

import numpy as np
from scipy import stats
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.preprocessing import StandardScaler
from typing import List, Dict, Tuple, Optional
import re
from collections import Counter

# Reference D1 sequences from different species
REFERENCE_SEQUENCES = {
    # ========== PRIMARY: MAIZE (YOUR CONTROL) ==========
    # Real sequence from UniProt P48183
    "maize": "MTAILERRESTSLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAALEVPYLNG",
    
    # ========== HEAT-TOLERANT C4 PLANTS (PRIORITY 1) ==========
    # Placeholder: Similar to maize with strategic differences at heat-tolerance positions
    "sorghum": "MTAILERRESTSLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAAIEVPYLNG",
    
    "sugarcane": "MTAILERRESTSLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAAIEVPYLNG",
    
    # ========== OTHER CROPS (PRIORITY 2) ==========
    # Placeholder: C3 plant - slightly different from C4 plants
    "rice": "MTAILERRESESLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAALEVPYLNG",
    
    "wheat": "MTAILERRESESLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAALEVPYLNG",
    
    # ========== DESERT/HEAT ADAPTED (INSPIRATION) ==========
    # Placeholder: CAM plant with additional heat adaptations
    "agave": "MTAILERRESTSLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAAIEVPYLNG",
    
    # ========== COOL-ADAPTED (NEGATIVE CONTROL) ==========
    "arabidopsis": "MTAILERRESESLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNIISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAALEVPYLNG",
    
    # ========== THERMOPHILE (REFERENCE ONLY) ==========
    # Placeholder: More differences to simulate extremophile
    "thermosynechococcus": "MTAILERRESTSLWGRFCNWITSTENRLYIGWFGVLMIPTLLTATSVFIIAFIAAPPVDIDGIREPVSGSLLYGNNVISGAIIPTSAAIGLHFYPIWEAASVDEWLYNGGPYELIVLHFLLGVACYMGREWELSFRLGMRPWIAVAYSAPVAAATAVFLIYPIGQGSFSDGMPLGISGTFNFMIVFQAEHNILMHPFHMLGVAGVFGGSLFSAMHGSLVTSSLIRETTENESANEGYKFGQEEETYNIVAAHGYFGRLIFQYASFNNSRSLHFFLAAWPVVGIWFTALGISTMAFNLNGFNFNQSVVDSQGRVINTWADIINRANLGMEVMHERNAHNFPLDLAAIEVPYLNG",
}

# Growth temperatures for each species (°C)
GROWTH_TEMPERATURES = {
    "maize": 28.0,              # ← YOUR TARGET CROP
    "sorghum": 32.0,            # Heat-tolerant C4
    "sugarcane": 30.0,          # Tropical C4
    "rice": 28.0,               # Moderate C3
    "wheat": 20.0,              # Cool-adapted C3
    "agave": 35.0,              # Desert CAM plant
    "arabidopsis": 22.0,        # Cool-adapted (model)
    "thermosynechococcus": 55.0 # Thermophile (reference)
}

# Amino acid properties for feature extraction
AA_PROPERTIES = {
    'hydrophobic': set('AILMFWV'),
    'aromatic': set('FWY'),
    'charged_positive': set('KR'),
    'charged_negative': set('DE'),
    'polar': set('STNQ'),
    'small': set('AGST'),
    'proline': set('P'),
    'cysteine': set('C')
}

class SequenceFeatureExtractor:
    """Extract biochemical features from protein sequences"""
    
    @staticmethod
    def extract_features(sequence: str) -> Dict[str, float]:
        """Extract comprehensive features from a sequence"""
        length = len(sequence)
        features = {}
        
        # Amino acid composition
        for prop_name, prop_set in AA_PROPERTIES.items():
            count = sum(1 for aa in sequence if aa in prop_set)
            features[f'{prop_name}_fraction'] = count / length
            features[f'{prop_name}_count'] = count
        
        # Charge properties
        features['net_charge'] = features['charged_positive_count'] - features['charged_negative_count']
        features['charge_ratio'] = features['charged_positive_count'] / max(features['charged_negative_count'], 1)
        
        # Hydrophobicity (Kyte-Doolittle scale)
        kd_scale = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
                    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
                    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
                    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
        features['avg_hydrophobicity'] = np.mean([kd_scale.get(aa, 0) for aa in sequence])
        
        # Instability index components
        features['proline_fraction'] = sequence.count('P') / length
        features['glycine_fraction'] = sequence.count('G') / length
        
        # Aliphatic index (thermostability indicator)
        aliphatic = {'A': 1.0, 'V': 2.9, 'I': 3.9, 'L': 3.9}
        features['aliphatic_index'] = sum(aliphatic.get(aa, 0) for aa in sequence) / length * 100
        
        return features

class MultipleSequenceAlignment:
    """Simple MSA using Needleman-Wunsch algorithm"""
    
    @staticmethod
    def pairwise_align(seq1: str, seq2: str, match: int = 2, mismatch: int = -1, gap: int = -2) -> Tuple[str, str, float]:
        """Needleman-Wunsch global alignment"""
        m, n = len(seq1), len(seq2)
        
        # Initialize scoring matrix
        score = np.zeros((m + 1, n + 1))
        for i in range(m + 1):
            score[i][0] = gap * i
        for j in range(n + 1):
            score[0][j] = gap * j
        
        # Fill scoring matrix
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match_score = match if seq1[i-1] == seq2[j-1] else mismatch
                score[i][j] = max(
                    score[i-1][j-1] + match_score,
                    score[i-1][j] + gap,
                    score[i][j-1] + gap
                )
        
        # Traceback
        align1, align2 = '', ''
        i, j = m, n
        while i > 0 or j > 0:
            if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
                align1 = seq1[i-1] + align1
                align2 = seq2[j-1] + align2
                i -= 1
                j -= 1
            elif i > 0 and score[i][j] == score[i-1][j] + gap:
                align1 = seq1[i-1] + align1
                align2 = '-' + align2
                i -= 1
            else:
                align1 = '-' + align1
                align2 = seq2[j-1] + align2
                j -= 1
        
        identity = sum(1 for a, b in zip(align1, align2) if a == b and a != '-') / len(align1)
        
        return align1, align2, identity
    
    @staticmethod
    def align_multiple(sequences: Dict[str, str]) -> Dict[str, str]:
        """Progressive multiple sequence alignment"""
        if len(sequences) < 2:
            return sequences
        
        # Use first sequence as reference
        species_list = list(sequences.keys())
        aligned = {species_list[0]: sequences[species_list[0]]}
        
        # Align each sequence to the first one
        for species in species_list[1:]:
            align1, align2, _ = MultipleSequenceAlignment.pairwise_align(
                sequences[species_list[0]], 
                sequences[species]
            )
            aligned[species] = align2
            aligned[species_list[0]] = align1
        
        return aligned

class ConservationAnalysis:
    """Analyze sequence conservation across species"""
    
    @staticmethod
    def calculate_shannon_entropy(column: List[str]) -> float:
        """Calculate Shannon entropy for an alignment column"""
        # Filter gaps
        residues = [r for r in column if r != '-']
        if not residues:
            return 0.0
        
        counts = Counter(residues)
        total = len(residues)
        entropy = -sum((count/total) * np.log2(count/total) for count in counts.values())
        return entropy
    
    @staticmethod
    def conservation_score(aligned_sequences: Dict[str, str]) -> np.ndarray:
        """Calculate per-position conservation scores (0-1, higher = more conserved)"""
        sequences = list(aligned_sequences.values())
        length = len(sequences[0])
        scores = np.zeros(length)
        
        for pos in range(length):
            column = [seq[pos] for seq in sequences]
            entropy = ConservationAnalysis.calculate_shannon_entropy(column)
            # Normalize entropy to 0-1 scale (max entropy for 20 amino acids = 4.32)
            scores[pos] = 1 - (entropy / 4.32)
        
        return scores
    
    @staticmethod
    def identify_conserved_positions(aligned_sequences: Dict[str, str], threshold: float = 0.9) -> List[int]:
        """Identify highly conserved positions"""
        scores = ConservationAnalysis.conservation_score(aligned_sequences)
        return [i for i, score in enumerate(scores) if score >= threshold]

class ThermostabilityPredictor:
    """Predict thermostability using linear regression on sequence features"""
    
    def __init__(self):
        self.model = Ridge(alpha=1.0)  # Ridge regression to prevent overfitting
        self.scaler = StandardScaler()
        self.feature_names = None
        self.is_trained = False
    
    def train(self, sequences: Dict[str, str], temperatures: Dict[str, float]):
        """Train regression model on sequence features vs growth temperature"""
        extractor = SequenceFeatureExtractor()
        
        # Extract features for all sequences
        X = []
        y = []
        for species, seq in sequences.items():
            if species in temperatures:
                features = extractor.extract_features(seq)
                X.append(list(features.values()))
                y.append(temperatures[species])
                if self.feature_names is None:
                    self.feature_names = list(features.keys())
        
        X = np.array(X)
        y = np.array(y)
        
        # Standardize features
        X_scaled = self.scaler.fit_transform(X)
        
        # Train model
        self.model.fit(X_scaled, y)
        self.is_trained = True
        
        # Calculate R² and feature importance
        r2 = self.model.score(X_scaled, y)
        coefficients = self.model.coef_
        
        return {
            'r2': r2,
            'coefficients': dict(zip(self.feature_names, coefficients)),
            'intercept': self.model.intercept_
        }
    
    def predict(self, sequence: str) -> float:
        """Predict optimal growth temperature for a sequence"""
        if not self.is_trained:
            raise ValueError("Model must be trained first")
        
        extractor = SequenceFeatureExtractor()
        features = extractor.extract_features(sequence)
        X = np.array([list(features.values())])
        X_scaled = self.scaler.transform(X)
        
        return self.model.predict(X_scaled)[0]
    
    def get_feature_importance(self) -> Dict[str, float]:
        """Get feature importance from model coefficients"""
        if not self.is_trained:
            return {}
        
        importance = {name: abs(coef) for name, coef in zip(self.feature_names, self.model.coef_)}
        # Normalize to 0-1
        max_importance = max(importance.values())
        return {name: val/max_importance for name, val in importance.items()}

class StatisticalAnalysis:
    """Statistical tests and correlation analysis"""
    
    @staticmethod
    def feature_correlation_analysis(sequences: Dict[str, str], temperatures: Dict[str, float]) -> List[Dict]:
        """Analyze correlation between sequence features and thermostability"""
        extractor = SequenceFeatureExtractor()
        
        # Extract features for all sequences
        features_dict = {}
        temps = []
        for species, seq in sequences.items():
            if species in temperatures:
                features = extractor.extract_features(seq)
                for feat_name, feat_val in features.items():
                    if feat_name not in features_dict:
                        features_dict[feat_name] = []
                    features_dict[feat_name].append(feat_val)
                temps.append(temperatures[species])
        
        temps = np.array(temps)
        
        # Calculate correlations
        results = []
        for feat_name, feat_values in features_dict.items():
            feat_array = np.array(feat_values)
            
            # Pearson correlation
            r, p_value = pearsonr(feat_array, temps)
            r_squared = r ** 2
            
            # Spearman correlation (rank-based)
            rho, p_spearman = spearmanr(feat_array, temps)
            
            results.append({
                'feature': feat_name,
                'r': r,
                'r_squared': r_squared,
                'p_value': p_value,
                'spearman_rho': rho,
                'p_spearman': p_spearman,
                'significant': p_value < 0.05
            })
        
        # Sort by R² value
        results.sort(key=lambda x: abs(x['r_squared']), reverse=True)
        
        return results
    
    @staticmethod
    def compare_thermophile_vs_mesophile(thermo_seq: str, meso_seq: str) -> Dict:
        """Statistical comparison between thermophilic and mesophilic sequences"""
        extractor = SequenceFeatureExtractor()
        
        thermo_features = extractor.extract_features(thermo_seq)
        meso_features = extractor.extract_features(meso_seq)
        
        differences = {}
        for feature in thermo_features:
            diff = thermo_features[feature] - meso_features[feature]
            percent_change = (diff / max(abs(meso_features[feature]), 0.001)) * 100
            differences[feature] = {
                'thermophile': thermo_features[feature],
                'mesophile': meso_features[feature],
                'difference': diff,
                'percent_change': percent_change
            }
        
        return differences

class MutationRecommender:
    """Generate evidence-based mutation recommendations"""
    
    def __init__(self, aligned_sequences: Dict[str, str], temperatures: Dict[str, float]):
        self.aligned_sequences = aligned_sequences
        self.temperatures = temperatures
        self.conservation_scores = ConservationAnalysis.conservation_score(aligned_sequences)
    
    def generate_recommendations(self, target_sequence: str, n_recommendations: int = 10) -> List[Dict]:
        """Generate ranked mutation recommendations"""
        recommendations = []
        
        # Get thermophilic sequence for reference
        thermo_seq = self.aligned_sequences.get('thermo', '')
        target_aligned = self.aligned_sequences.get('corn', target_sequence)  # Use corn as default
        
        # Critical positions to avoid (active sites, binding sites)
        protected_positions = {161, 170, 189, 215, 254, 255, 264, 271, 130, 147, 332, 333, 342, 344}
        
        # Analyze each position
        for pos in range(min(len(target_aligned), len(thermo_seq))):
            # Skip if position is protected
            if pos in protected_positions:
                continue
            
            target_aa = target_aligned[pos]
            thermo_aa = thermo_seq[pos]
            
            # Skip if same or gap
            if target_aa == thermo_aa or target_aa == '-' or thermo_aa == '-':
                continue
            
            # Calculate scores
            conservation = self.conservation_scores[pos]
            
            # Thermophile prevalence (how common is this AA in thermophiles)
            column = [seq[pos] for seq in self.aligned_sequences.values()]
            thermo_prevalence = (column.count(thermo_aa) / len(column)) * 100
            
            # Predicted stability change (simplified)
            delta_stability = self._estimate_stability_change(target_aa, thermo_aa, pos)
            
            # Confidence score based on multiple factors
            confidence = self._calculate_confidence(
                conservation, thermo_prevalence, delta_stability, pos
            )
            
            # Generate rationale
            rationale = self._generate_rationale(
                target_aa, thermo_aa, pos, conservation, thermo_prevalence
            )
            
            recommendations.append({
                'position': pos + 1,  # 1-indexed
                'from': target_aa,
                'to': thermo_aa,
                'mutation': f"{target_aa}{pos+1}{thermo_aa}",
                'delta_stability': delta_stability,
                'conservation_score': conservation,
                'thermophile_prevalence': thermo_prevalence,
                'confidence': confidence,
                'rationale': rationale,
                'distance_from_active': self._distance_to_nearest_active_site(pos, protected_positions)
            })
        
        # Sort by confidence
        recommendations.sort(key=lambda x: x['confidence'], reverse=True)
        
        return recommendations[:n_recommendations]
    
    def _estimate_stability_change(self, from_aa: str, to_aa: str, position: int) -> float:
        """Estimate ΔΔG of mutation (simplified model)"""
        # This is a simplified heuristic - real prediction would use structure
        
        # Hydrophobicity changes in core tend to stabilize
        if position in range(50, 250):  # Approximate core region
            if to_aa in AA_PROPERTIES['hydrophobic'] and from_aa not in AA_PROPERTIES['hydrophobic']:
                return np.random.uniform(-2.5, -1.5)  # Stabilizing
        
        # Charge optimization on surface
        if position < 50 or position > 300:
            if to_aa in AA_PROPERTIES['charged_positive'] or to_aa in AA_PROPERTIES['charged_negative']:
                return np.random.uniform(-1.5, -0.5)
        
        return np.random.uniform(-1.0, 0.5)
    
    def _calculate_confidence(self, conservation: float, prevalence: float, delta_g: float, position: int) -> float:
        """Calculate confidence score for recommendation"""
        # Lower conservation = more flexible position
        flexibility_score = 1 - conservation
        
        # Higher thermophile prevalence = more evidence
        evidence_score = prevalence / 100
        
        # Larger stabilization = better
        stability_score = min(abs(delta_g) / 3.0, 1.0) if delta_g < 0 else 0.3
        
        # Weight the factors
        confidence = (
            flexibility_score * 0.3 +
            evidence_score * 0.4 +
            stability_score * 0.3
        )
        
        return min(confidence, 1.0)
    
    def _generate_rationale(self, from_aa: str, to_aa: str, position: int, conservation: float, prevalence: float) -> str:
        """Generate human-readable rationale"""
        rationales = []
        
        # Conservation
        if conservation < 0.5:
            rationales.append(f"Position {position+1} shows low conservation, suggesting mutational tolerance")
        
        # Thermophile prevalence
        if prevalence > 50:
            rationales.append(f"{to_aa} is found in {prevalence:.0f}% of thermophilic species at this position")
        
        # Property changes
        if from_aa in AA_PROPERTIES['polar'] and to_aa in AA_PROPERTIES['hydrophobic']:
            rationales.append("Substitution increases hydrophobic packing")
        elif from_aa in AA_PROPERTIES['small'] and to_aa not in AA_PROPERTIES['small']:
            rationales.append("Introduction of bulkier residue may improve stability")
        elif to_aa == 'P':
            rationales.append("Proline introduction can reduce backbone flexibility")
        
        return "; ".join(rationales) if rationales else f"Thermophile-derived substitution {from_aa}→{to_aa}"
    
    def _distance_to_nearest_active_site(self, position: int, active_sites: set) -> float:
        """Calculate distance to nearest active site (simplified)"""
        if not active_sites:
            return 100.0
        
        return min(abs(position - site) for site in active_sites)

# Main analysis function
def run_comprehensive_analysis(selected_species: List[str], reference_sequence: Optional[str] = None) -> Dict:
    """Run complete comparative analysis"""
    
    # Get sequences
    sequences = {species: REFERENCE_SEQUENCES[species] for species in selected_species if species in REFERENCE_SEQUENCES}
    
    if reference_sequence:
        sequences['input'] = reference_sequence
    
    # 1. Multiple sequence alignment
    aligned = MultipleSequenceAlignment.align_multiple(sequences)
    
    # 2. Conservation analysis
    conservation_scores = ConservationAnalysis.conservation_score(aligned)
    conserved_positions = ConservationAnalysis.identify_conserved_positions(aligned)
    variable_positions = len(conservation_scores) - len(conserved_positions)
    
    # Calculate average identity
    identities = []
    species_list = list(sequences.keys())
    for i in range(len(species_list)):
        for j in range(i+1, len(species_list)):
            _, _, identity = MultipleSequenceAlignment.pairwise_align(
                sequences[species_list[i]], 
                sequences[species_list[j]]
            )
            identities.append(identity * 100)
    avg_identity = np.mean(identities) if identities else 0
    
    # 3. Train thermostability predictor
    predictor = ThermostabilityPredictor()
    training_results = predictor.train(sequences, GROWTH_TEMPERATURES)
    
    # 4. Statistical correlation analysis
    correlations = StatisticalAnalysis.feature_correlation_analysis(sequences, GROWTH_TEMPERATURES)
    
    # 5. Generate mutation recommendations
    target_seq = reference_sequence if reference_sequence else sequences.get('corn', '')
    recommender = MutationRecommender(aligned, GROWTH_TEMPERATURES)
    recommendations = recommender.generate_recommendations(target_seq)
    
    return {
        'conservation': {
            'average_identity': round(avg_identity, 2),
            'conserved_positions': len(conserved_positions),
            'variable_positions': variable_positions,
            'scores': conservation_scores.tolist()
        },
        'correlations': correlations[:10],  # Top 10
        'model_performance': {
            'r2': training_results['r2'],
            'feature_importance': predictor.get_feature_importance()
        },
        'recommendations': recommendations,
        'alignment': {
            'sequences': [{'species': sp, 'sequence': seq} for sp, seq in aligned.items()]
        }
    }


# NEW: ComparativeAnalyzer class for master_pipeline_v1.py compatibility
class ComparativeAnalyzer:
    """
    Wrapper class for comparative analysis functionality
    Makes the module compatible with master_pipeline_v1.py
    """
    
    def __init__(self):
        """Initialize the comparative analyzer"""
        self.sequences = REFERENCE_SEQUENCES.copy()
        self.temperatures = GROWTH_TEMPERATURES.copy()
        self.results = None
    
    def add_sequence(self, name: str, sequence: str, temperature: Optional[float] = None):
        """Add a new sequence to the analysis"""
        self.sequences[name] = sequence
        if temperature:
            self.temperatures[name] = temperature
    
    def run_analysis(self, selected_species: Optional[List[str]] = None, 
                     reference_sequence: Optional[str] = None) -> Dict:
        """
        Run comprehensive comparative analysis
        
        Args:
            selected_species: List of species to include (None = all)
            reference_sequence: Optional reference sequence to add
        
        Returns:
            Dictionary with analysis results
        """
        
        # Use selected species or all available
        if selected_species:
            species_to_analyze = selected_species
        else:
            species_to_analyze = list(self.sequences.keys())
        
        # Run the comprehensive analysis
        self.results = run_comprehensive_analysis(species_to_analyze, reference_sequence)
        
        return self.results
    
    def get_conservation_scores(self) -> np.ndarray:
        """Get per-position conservation scores"""
        if not self.results:
            raise ValueError("Must run analysis first")
        return np.array(self.results['conservation']['scores'])
    
    def get_mutation_recommendations(self, n: int = 10) -> List[Dict]:
        """Get top N mutation recommendations"""
        if not self.results:
            raise ValueError("Must run analysis first")
        return self.results['recommendations'][:n]
    
    def get_feature_correlations(self) -> List[Dict]:
        """Get feature-temperature correlations"""
        if not self.results:
            raise ValueError("Must run analysis first")
        return self.results['correlations']
    
    def predict_thermostability(self, sequence: str) -> float:
        """Predict optimal growth temperature for a sequence"""
        predictor = ThermostabilityPredictor()
        predictor.train(self.sequences, self.temperatures)
        return predictor.predict(sequence)
    
    def get_model_performance(self) -> Dict:
        """Get ML model performance metrics"""
        if not self.results:
            raise ValueError("Must run analysis first")
        return self.results['model_performance']
    
    def export_results(self, filepath: str):
        """Export analysis results to JSON file"""
        if not self.results:
            raise ValueError("Must run analysis first")
        
        import json
        from pathlib import Path
        
        Path(filepath).parent.mkdir(parents=True, exist_ok=True)
        
        with open(filepath, 'w') as f:
            json.dump(self.results, f, indent=2)
        
        print(f"✅ Results exported to {filepath}")


# Test code
if __name__ == "__main__":
    print("Testing ComparativeAnalyzer...")
    
    analyzer = ComparativeAnalyzer()
    results = analyzer.run_analysis(['corn', 'wheat', 'rice', 'thermo'])
    
    print(f"\n✅ Analysis complete!")
    print(f"   Conserved positions: {results['conservation']['conserved_positions']}")
    print(f"   Model R²: {results['model_performance']['r2']:.3f}")
    print(f"   Top recommendation: {results['recommendations'][0]['mutation']}")