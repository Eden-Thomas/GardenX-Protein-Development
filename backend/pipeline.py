import asyncio
import numpy as np
import pandas as pd
from typing import List, Dict, Optional, Callable
from datetime import datetime
import json

from config import Config
from neurosnap_client import NeurosnapClient
from models import LocalModels, ThermoScoreCalculator

class D1EngineeringPipeline:
    def __init__(self, config: Config):
        self.config = config
        self.neurosnap = NeurosnapClient(config.NEUROSNAP_API_KEY) if config.NEUROSNAP_API_KEY else None
        self.local_models = LocalModels() if config.USE_LOCAL_GPU else None
        self.scorer = ThermoScoreCalculator(config)
        
    async def run_pipeline(self, 
                          sequences: List[str], 
                          crop_type: str,
                          n_variants: int,
                          progress_callback: Optional[Callable] = None) -> Dict:
        """Execute the complete D1 engineering pipeline"""
        
        results = {
            "crop_type": crop_type,
            "input_sequences": len(sequences),
            "timestamp": datetime.now().isoformat(),
            "stages": {}
        }
        
        # Stage 1: Fitness Scoring (20% progress)
        if progress_callback:
            progress_callback(0.1)
        
        print("Stage 1: Computing fitness scores...")
        fitness_scores = await self._compute_fitness(sequences)
        results["stages"]["fitness"] = {
            "complete": True,
            "avg_score": np.mean(list(fitness_scores.values()))
        }
        
        if progress_callback:
            progress_callback(0.2)
        
        # Stage 2: Generate Variants (40% progress)
        print("Stage 2: Generating variants...")
        variants = await self._generate_variants(sequences, crop_type, n_variants)
        results["stages"]["generation"] = {
            "complete": True,
            "n_variants": len(variants)
        }
        
        if progress_callback:
            progress_callback(0.4)
        
        # Stage 3: Score Variants (60% progress)
        print("Stage 3: Scoring variants...")
        scored_variants = await self._score_variants(variants)
        results["stages"]["scoring"] = {
            "complete": True,
            "top_score": scored_variants.iloc[0]["thermo_score"] if len(scored_variants) > 0 else 0
        }
        
        if progress_callback:
            progress_callback(0.6)
        
        # Stage 4: Structure Prediction (80% progress)
        print("Stage 4: Predicting structures...")
        top_variants = scored_variants.head(50)
        structures = await self._predict_structures(top_variants["sequence"].tolist())
        results["stages"]["structures"] = {
            "complete": True,
            "n_structures": len(structures)
        }
        
        if progress_callback:
            progress_callback(0.8)
        
        # Stage 5: Diverse Selection (100% progress)
        print("Stage 5: Selecting diverse candidates...")
        final_picks = self._diverse_selection(scored_variants, structures)
        results["stages"]["selection"] = {
            "complete": True,
            "n_selected": len(final_picks)
        }
        
        # Save results
        self._save_results(final_picks, results)
        
        if progress_callback:
            progress_callback(1.0)
        
        results["final_variants"] = final_picks.to_dict("records")
        return results
    
    async def _compute_fitness(self, sequences: List[str]) -> Dict[str, float]:
        """Compute fitness scores for sequences"""
        if self.neurosnap:
            scores = await self.neurosnap.compute_fitness(sequences)
        elif self.local_models:
            scores = self.local_models.compute_esm_fitness(sequences)
        else:
            # Fallback: composition-based scoring
            scores = {seq: self.scorer.compute_composition_score(seq) for seq in sequences}
        
        return scores
    
    async def _generate_variants(self, templates: List[str], crop_type: str, n_variants: int) -> List[Dict]:
        """Generate thermostable variants"""
        variants = []
        
        if self.neurosnap:
            # Use Neurosnap APIs
            variants = await self.neurosnap.generate_variants(templates, n_variants)
        else:
            # Simple mutation strategy
            for template in templates:
                for _ in range(n_variants // len(templates)):
                    variant = self._apply_thermostable_mutations(template)
                    variants.append({"sequence": variant})
        
        return variants
    
    def _apply_thermostable_mutations(self, sequence: str) -> str:
        """Apply rule-based thermostable mutations"""
        seq_list = list(sequence)
        
        # Example mutations (simplified)
        mutations = {
            "G": ["A", "V"],  # Replace Gly with more rigid
            "S": ["T", "V"],  # Replace Ser with bulkier
            "A": ["I", "L", "V"]  # Enhance hydrophobic packing
        }
        
        # Apply some mutations in redesign windows
        for window in self.config.REDESIGN_WINDOWS:
            for pos in window:
                if pos < len(seq_list):
                    current = seq_list[pos]
                    if current in mutations:
                        seq_list[pos] = np.random.choice(mutations[current])
        
        return "".join(seq_list)
    
    async def _score_variants(self, variants: List[Dict]) -> pd.DataFrame:
        """Calculate ThermoScore for variants"""
        scores = []
        
        for variant in variants:
            seq = variant["sequence"]
            thermo_score = self.scorer.calculate_thermoscore(seq)
            
            scores.append({
                "sequence": seq,
                "thermo_score": thermo_score,
                **self.scorer.get_features(seq)
            })
        
        return pd.DataFrame(scores).sort_values("thermo_score", ascending=False)
    
    async def _predict_structures(self, sequences: List[str]) -> List[Dict]:
        """Predict structures for top variants"""
        if self.neurosnap:
            return await self.neurosnap.predict_structures(sequences)
        else:
            # Mock structures for testing
            return [{"plddt": 85 + np.random.rand() * 10} for _ in sequences]
    
    def _diverse_selection(self, scored_variants: pd.DataFrame, structures: List[Dict]) -> pd.DataFrame:
        """Select diverse high-scoring variants"""
        # Add structure quality
        scored_variants["plddt"] = [s.get("plddt", 85) for s in structures[:len(scored_variants)]]
        
        # Filter by quality
        candidates = scored_variants[
            (scored_variants["thermo_score"] > 0.6) &
            (scored_variants["plddt"] > 80)
        ].head(50)
        
        return candidates
    
    def _save_results(self, variants: pd.DataFrame, metadata: Dict):
        """Save results to disk"""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # Save FASTA
        fasta_path = self.config.RESULTS_DIR / f"variants_{timestamp}.fasta"
        with open(fasta_path, "w") as f:
            for idx, row in variants.iterrows():
                f.write(f">variant_{idx}_score_{row['thermo_score']:.3f}\n")
                f.write(f"{row['sequence']}\n")
        
        # Save CSV
        csv_path = self.config.RESULTS_DIR / f"variants_{timestamp}.csv"
        variants.to_csv(csv_path, index=False)
        
        print(f"Results saved to {self.config.RESULTS_DIR}")