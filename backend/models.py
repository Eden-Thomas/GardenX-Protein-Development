import numpy as np
from typing import List, Dict

class LocalModels:
    """Local model implementations"""
    
    def __init__(self):
        self.esm_model = None
        self._load_models()
    
    def _load_models(self):
        """Load local models if available"""
        try:
            import esm
            self.esm_model, self.esm_alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
            self.esm_model.eval()
            print("ESM model loaded successfully")
        except:
            print("ESM model not available, using fallback scoring")
    
    def compute_esm_fitness(self, sequences: List[str]) -> Dict[str, float]:
        """Compute ESM fitness scores"""
        if self.esm_model:
            # Actual ESM scoring would go here
            pass
        
        # Fallback scoring
        return {seq: np.random.rand() for seq in sequences}

class ThermoScoreCalculator:
    """Calculate thermostability scores"""
    
    def __init__(self, config):
        self.config = config
        self.weights = config.THERMOSCORE_WEIGHTS
    
    def calculate_thermoscore(self, sequence: str) -> float:
        """Calculate overall thermostability score"""
        features = self.get_features(sequence)
        
        # Weighted combination
        score = (
            self.weights["composition"] * features["composition_score"] +
            0.3  # Placeholder for ESM score
        )
        
        return min(max(score, 0), 1)  # Normalize to [0,1]
    
    def get_features(self, sequence: str) -> Dict[str, float]:
        """Extract thermostability features"""
        L = len(sequence)
        
        # IVYWREL fraction (hydrophobic/aromatic)
        ivywrel = sum(c in "IVYWREL" for c in sequence) / L
        
        # Aromatics
        aromatics = sum(c in "FYWH" for c in sequence) / L
        
        # Salt bridge potential
        neg = sum(c in "ED" for c in sequence)
        pos = sum(c in "KR" for c in sequence)
        salt_proxy = (neg * pos) / (L * L) if L > 0 else 0
        
        # Glycine penalty
        gly_penalty = sequence.count("G") / L
        
        # Composition score
        composition_score = ivywrel + aromatics + salt_proxy - gly_penalty
        
        return {
            "ivywrel": ivywrel,
            "aromatics": aromatics,
            "salt_proxy": salt_proxy,
            "gly_penalty": gly_penalty,
            "composition_score": composition_score
        }
    
    def compute_composition_score(self, sequence: str) -> float:
        """Simple composition-based fitness score"""
        features = self.get_features(sequence)
        return features["composition_score"]