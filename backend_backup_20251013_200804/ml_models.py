"""
ml_models.py
Machine learning models for predicting mutation effects on protein stability
"""

import numpy as np
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score
import joblib
from typing import Dict, List, Tuple, Optional

class MutationEffectPredictor:
    """Ensemble model for predicting ΔΔG of mutations"""
    
    def __init__(self):
        self.models = {
            'random_forest': RandomForestRegressor(n_estimators=100, random_state=42),
            'gradient_boosting': GradientBoostingRegressor(n_estimators=100, random_state=42),
            'neural_network': MLPRegressor(hidden_layer_sizes=(100, 50), max_iter=1000, random_state=42)
        }
        self.scaler = StandardScaler()
        self.is_trained = False
        self.ensemble_weights = None
    
    def encode_mutation(self, from_aa: str, to_aa: str, position: int, sequence_length: int) -> np.ndarray:
        """Encode mutation features for ML model"""
        features = []
        
        # One-hot encoding for amino acids (20 dimensions each)
        aa_order = 'ACDEFGHIKLMNPQRSTVWY'
        from_encoding = [1 if aa == from_aa else 0 for aa in aa_order]
        to_encoding = [1 if aa == to_aa else 0 for aa in aa_order]
        features.extend(from_encoding)
        features.extend(to_encoding)
        
        # Position features
        features.append(position / sequence_length)  # Relative position
        features.append(1 if position < sequence_length * 0.2 else 0)  # N-terminal
        features.append(1 if position > sequence_length * 0.8 else 0)  # C-terminal
        features.append(1 if 0.3 < position/sequence_length < 0.7 else 0)  # Core region
        
        # Physicochemical property changes
        hydrophobic = set('AILMFWV')
        aromatic = set('FWY')
        charged_pos = set('KR')
        charged_neg = set('DE')
        polar = set('STNQ')
        
        features.append(1 if (from_aa in hydrophobic) != (to_aa in hydrophobic) else 0)
        features.append(1 if (from_aa in aromatic) != (to_aa in aromatic) else 0)
        features.append(1 if (from_aa in charged_pos) != (to_aa in charged_pos) else 0)
        features.append(1 if (from_aa in charged_neg) != (to_aa in charged_neg) else 0)
        features.append(1 if (from_aa in polar) != (to_aa in polar) else 0)
        
        # Size change
        size_scale = {'G': 1, 'A': 2, 'S': 2, 'C': 2, 'D': 3, 'P': 3, 'N': 3, 'T': 3,
                      'E': 4, 'V': 4, 'Q': 4, 'H': 4, 'M': 4, 'I': 4, 'L': 4, 'K': 4,
                      'R': 5, 'F': 5, 'Y': 5, 'W': 6}
        size_change = size_scale.get(to_aa, 3) - size_scale.get(from_aa, 3)
        features.append(size_change / 5)  # Normalize
        
        # Proline effects
        features.append(1 if to_aa == 'P' and from_aa != 'P' else 0)
        features.append(1 if from_aa == 'P' and to_aa != 'P' else 0)
        
        return np.array(features)
    
    def train_with_synthetic_data(self, n_samples: int = 1000):
        """Train models using synthetic data based on known mutation principles"""
        # Generate synthetic training data
        X = []
        y = []
        
        aa_list = list('ACDEFGHIKLMNPQRSTVWY')
        sequence_length = 350  # Approximate D1 length
        
        for _ in range(n_samples):
            pos = np.random.randint(10, sequence_length - 10)
            from_aa = np.random.choice(aa_list)
            to_aa = np.random.choice([a for a in aa_list if a != from_aa])
            
            # Encode features
            features = self.encode_mutation(from_aa, to_aa, pos, sequence_length)
            X.append(features)
            
            # Simulate ΔΔG based on rules
            ddg = self._simulate_ddg(from_aa, to_aa, pos, sequence_length)
            y.append(ddg)
        
        X = np.array(X)
        y = np.array(y)
        
        # Scale features
        X_scaled = self.scaler.fit_transform(X)
        
        # Train each model and calculate weights based on cross-validation
        cv_scores = {}
        for name, model in self.models.items():
            model.fit(X_scaled, y)
            scores = cross_val_score(model, X_scaled, y, cv=5, scoring='r2')
            cv_scores[name] = np.mean(scores)
        
        # Calculate ensemble weights (proportional to performance)
        total_score = sum(max(0, s) for s in cv_scores.values())
        self.ensemble_weights = {
            name: max(0, score) / total_score 
            for name, score in cv_scores.items()
        }
        
        self.is_trained = True
        
        return {
            'cv_scores': cv_scores,
            'ensemble_weights': self.ensemble_weights
        }
    
    def _simulate_ddg(self, from_aa: str, to_aa: str, position: int, seq_length: int) -> float:
        """Simulate ΔΔG based on biophysical principles"""
        ddg = 0.0
        
        # Core region mutations
        if 0.3 < position/seq_length < 0.7:
            # Hydrophobic in core is stabilizing
            if to_aa in 'ILMFV' and from_aa not in 'ILMFV':
                ddg -= np.random.uniform(1.5, 3.0)
            # Polar in core is destabilizing
            elif to_aa in 'STNQ' and from_aa not in 'STNQ':
                ddg += np.random.uniform(1.0, 2.5)
        
        # Surface mutations
        else:
            # Charged residues on surface can stabilize
            if to_aa in 'KRED' and from_aa not in 'KRED':
                ddg -= np.random.uniform(0.5, 1.5)
        
        # Proline effects
        if to_aa == 'P' and from_aa != 'P':
            # Proline reduces backbone flexibility
            ddg -= np.random.uniform(0.5, 1.5)
        elif from_aa == 'P' and to_aa != 'P':
            ddg += np.random.uniform(0.5, 1.5)
        
        # Aromatic packing
        if to_aa in 'FWY' and from_aa in 'FWY':
            ddg -= np.random.uniform(0.3, 0.8)
        
        # Size mismatches
        size_scale = {'G': 1, 'A': 2, 'S': 2, 'V': 3, 'I': 4, 'L': 4, 'W': 6}
        from_size = size_scale.get(from_aa, 3)
        to_size = size_scale.get(to_aa, 3)
        if abs(to_size - from_size) > 2:
            ddg += np.random.uniform(0.5, 1.5)
        
        # Add noise
        ddg += np.random.normal(0, 0.3)
        
        return ddg
    
    def predict(self, from_aa: str, to_aa: str, position: int, sequence_length: int = 350) -> Dict:
        """Predict ΔΔG using ensemble of models"""
        if not self.is_trained:
            self.train_with_synthetic_data()
        
        features = self.encode_mutation(from_aa, to_aa, position, sequence_length)
        X = self.scaler.transform([features])
        
        # Get predictions from all models
        predictions = {}
        for name, model in self.models.items():
            predictions[name] = model.predict(X)[0]
        
        # Ensemble prediction (weighted average)
        ensemble_pred = sum(
            pred * self.ensemble_weights[name] 
            for name, pred in predictions.items()
        )
        
        # Estimate uncertainty (std of predictions)
        uncertainty = np.std(list(predictions.values()))
        
        return {
            'predicted_ddg': ensemble_pred,
            'uncertainty': uncertainty,
            'individual_predictions': predictions,
            'confidence': self._calculate_confidence(ensemble_pred, uncertainty)
        }
    
    def _calculate_confidence(self, ddg: float, uncertainty: float) -> float:
        """Calculate confidence score based on prediction magnitude and uncertainty"""
        # Higher magnitude = more confident
        magnitude_score = min(abs(ddg) / 3.0, 1.0)
        
        # Lower uncertainty = more confident
        uncertainty_score = max(0, 1 - uncertainty / 2.0)
        
        return (magnitude_score + uncertainty_score) / 2
    
    def predict_multiple(self, mutations: List[Tuple[str, str, int]], sequence_length: int = 350) -> List[Dict]:
        """Predict effects for multiple mutations"""
        results = []
        for from_aa, to_aa, position in mutations:
            pred = self.predict(from_aa, to_aa, position, sequence_length)
            pred['mutation'] = f"{from_aa}{position}{to_aa}"
            results.append(pred)
        
        return results

class VariantScorer:
    """Score complete variant sequences for thermostability"""
    
    def __init__(self):
        self.mutation_predictor = MutationEffectPredictor()
        self.mutation_predictor.train_with_synthetic_data()
    
    def score_variant(self, reference_seq: str, variant_seq: str) -> Dict:
        """Score a variant against reference sequence"""
        if len(reference_seq) != len(variant_seq):
            raise ValueError("Sequences must be same length")
        
        # Identify all mutations
        mutations = []
        for i, (ref_aa, var_aa) in enumerate(zip(reference_seq, variant_seq)):
            if ref_aa != var_aa:
                mutations.append((ref_aa, var_aa, i + 1))
        
        # Predict effect of each mutation
        mutation_effects = self.mutation_predictor.predict_multiple(mutations, len(reference_seq))
        
        # Calculate cumulative ΔΔG
        total_ddg = sum(m['predicted_ddg'] for m in mutation_effects)
        avg_confidence = np.mean([m['confidence'] for m in mutation_effects])
        
        # Estimate ΔTm (rough conversion: -1 kcal/mol ≈ +2°C)
        estimated_delta_tm = -total_ddg * 2.0
        
        return {
            'num_mutations': len(mutations),
            'mutations': mutations,
            'mutation_details': mutation_effects,
            'total_ddg': total_ddg,
            'estimated_delta_tm': estimated_delta_tm,
            'avg_confidence': avg_confidence,
            'thermoscore': self._calculate_thermoscore(total_ddg, avg_confidence)
        }
    
    def _calculate_thermoscore(self, ddg: float, confidence: float) -> float:
        """Calculate overall thermostability score (0-1)"""
        # Normalize ΔΔG to 0-1 scale (assume max improvement ~10 kcal/mol)
        stability_component = min(max(-ddg / 10.0, 0), 1)
        
        # Weight by confidence
        score = stability_component * confidence
        
        return score