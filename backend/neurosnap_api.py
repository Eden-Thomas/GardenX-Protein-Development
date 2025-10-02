"""
neurosnap_api.py
Complete NeuroSnap API integration for D1 protein engineering
"""

import requests
import time
import json
from typing import Dict, List, Optional, Tuple
from requests_toolbelt.multipart.encoder import MultipartEncoder
import numpy as np
from scipy import stats

class NeuroSnapAPI:
    """NeuroSnap API client for protein engineering"""
    
    def __init__(self, api_key: str):
        self.api_key = api_key
        self.base_url = "https://neurosnap.ai/api"
        self.headers = {"X-API-KEY": api_key}
    
    def _submit_job(self, service_name: str, fields: Dict) -> str:
        """Submit a job to NeuroSnap"""
        multipart_data = MultipartEncoder(fields=fields)
        
        response = requests.post(
            f"{self.base_url}/job/submit/{service_name}",
            headers={
                **self.headers,
                "Content-Type": multipart_data.content_type
            },
            data=multipart_data
        )
        
        if response.status_code == 200:
            return response.json()  # Returns job ID
        else:
            raise Exception(f"Job submission failed: {response.text}")
    
    def _check_status(self, job_id: str) -> str:
        """Check job status"""
        response = requests.get(
            f"{self.base_url}/job/status/{job_id}",
            headers=self.headers
        )
        return response.json()
    
    def _get_results(self, job_id: str) -> Dict:
        """Get job results"""
        response = requests.get(
            f"{self.base_url}/job/output/{job_id}",
            headers=self.headers
        )
        return response.json()
    
    def _wait_for_completion(self, job_id: str, timeout: int = 3600) -> Dict:
        """Wait for job to complete and return results"""
        start_time = time.time()
        
        while time.time() - start_time < timeout:
            status = self._check_status(job_id)
            
            if status == "completed":
                return self._get_results(job_id)
            elif status == "failed":
                raise Exception(f"Job {job_id} failed")
            
            time.sleep(10)  # Check every 10 seconds
        
        raise TimeoutError(f"Job {job_id} timed out after {timeout} seconds")
    
    # ========== AlphaFold2 Integration ==========
    
    def predict_structure(self, 
                         sequence: str, 
                         job_name: str = "D1_structure",
                         use_templates: bool = False) -> Dict:
        """
        Predict protein structure using AlphaFold2
        
        Returns:
            - pdb_file: PDB structure
            - plddt_scores: Per-residue confidence scores
            - pae_matrix: Predicted aligned error matrix
        """
        fields = {
            "Sequence": sequence,
            "MSA Mode": "MMseqs2 (UniRef+Environmental)",
            "Use Templates": "true" if use_templates else "false",
            "Max Recycles": "3",
            "Model Type": "auto",
            "Relax": "false",
            "Note": job_name
        }
        
        job_id = self._submit_job("AlphaFold2", fields)
        results = self._wait_for_completion(job_id)
        
        return {
            'job_id': job_id,
            'pdb_structure': results.get('ranked_0.pdb'),
            'plddt_scores': self._parse_plddt(results),
            'pae_matrix': results.get('pae.json'),
            'confidence': np.mean(self._parse_plddt(results))
        }
    
    def _parse_plddt(self, results: Dict) -> np.ndarray:
        """Extract pLDDT scores from AlphaFold2 results"""
        # Parse PDB file to extract B-factor column (contains pLDDT)
        pdb_content = results.get('ranked_0.pdb', '')
        plddt_scores = []
        
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                # B-factor is in columns 61-66
                plddt = float(line[60:66].strip())
                plddt_scores.append(plddt)
        
        return np.array(plddt_scores)
    
    # ========== NeuroFold Integration ==========
    
    def optimize_thermostability(self,
                                 sequence: str,
                                 num_designs: int = 10,
                                 job_name: str = "D1_thermostability") -> List[Dict]:
        """
        Optimize protein for thermostability using NeuroFold
        
        Returns list of optimized variants with predictions
        """
        fields = {
            "Sequence": sequence,
            "Number of Designs": str(num_designs),
            "Optimization Target": "Thermostability",
            "Note": job_name
        }
        
        job_id = self._submit_job("NeuroFold", fields)
        results = self._wait_for_completion(job_id, timeout=7200)  # 2 hours
        
        # Parse NeuroFold results
        variants = []
        for i in range(num_designs):
            variant_key = f"design_{i}"
            if variant_key in results:
                variants.append({
                    'sequence': results[variant_key]['sequence'],
                    'predicted_tm_increase': results[variant_key].get('delta_tm', 0),
                    'confidence': results[variant_key].get('confidence', 0),
                    'mutations': self._identify_mutations(sequence, results[variant_key]['sequence'])
                })
        
        return variants
    
    def _identify_mutations(self, ref_seq: str, var_seq: str) -> List[str]:
        """Identify mutations between reference and variant"""
        mutations = []
        for i, (ref, var) in enumerate(zip(ref_seq, var_seq)):
            if ref != var:
                mutations.append(f"{ref}{i+1}{var}")
        return mutations
    
    # ========== ProteinMPNN Integration ==========
    
    def inverse_fold_sequence(self,
                              pdb_structure: str,
                              num_samples: int = 10,
                              temperature: float = 0.1) -> List[Dict]:
        """
        Generate sequence designs from structure using ProteinMPNN
        """
        fields = {
            "Input Structure": json.dumps([{"type": "pdb", "data": pdb_structure}]),
            "Number of Designs": str(num_samples),
            "Temperature": str(temperature),
            "Note": "D1_inverse_folding"
        }
        
        job_id = self._submit_job("ProteinMPNN", fields)
        results = self._wait_for_completion(job_id)
        
        designs = []
        for i in range(num_samples):
            design_key = f"design_{i}"
            if design_key in results:
                designs.append({
                    'sequence': results[design_key]['sequence'],
                    'score': results[design_key].get('score', 0),
                    'recovery_rate': results[design_key].get('recovery', 0)
                })
        
        return designs
    
    # ========== RFdiffusion Integration ==========
    
    def generate_stable_structure(self,
                                   target_length: int = 350,
                                   num_designs: int = 5) -> List[Dict]:
        """
        Generate de novo stable protein structures using RFdiffusion
        """
        fields = {
            "Design Type": "Unconditional",
            "Target Length": str(target_length),
            "Number of Designs": str(num_designs),
            "Note": "D1_structure_generation"
        }
        
        job_id = self._submit_job("RFdiffusion", fields)
        results = self._wait_for_completion(job_id, timeout=3600)
        
        structures = []
        for i in range(num_designs):
            structure_key = f"design_{i}.pdb"
            if structure_key in results:
                structures.append({
                    'pdb_structure': results[structure_key],
                    'design_id': i
                })
        
        return structures


class EnhancedStatisticalAnalyzer:
    """
    Enhanced statistical analysis combining local predictions with NeuroSnap validation
    """
    
    def __init__(self, neurosnap_api: NeuroSnapAPI):
        self.neurosnap = neurosnap_api
    
    def validate_predictions_with_structure(self,
                                           predictions: List[Dict],
                                           reference_sequence: str) -> Dict:
        """
        Validate mutation predictions using AlphaFold2 structure predictions
        
        Returns:
            - correlation_r2: R² between predicted ΔΔG and structural confidence changes
            - p_value: Statistical significance
            - validated_predictions: Predictions with structure-based confidence
        """
        print("Running AlphaFold2 validation...")
        
        # Predict reference structure
        ref_structure = self.neurosnap.predict_structure(
            reference_sequence,
            job_name="D1_reference_structure"
        )
        ref_plddt = ref_structure['plddt_scores']
        
        validated_predictions = []
        structure_confidence_changes = []
        predicted_ddgs = []
        
        # Validate top 10 mutations (to save API costs)
        for pred in predictions[:10]:
            mutation_pos = pred['position'] - 1  # Convert to 0-indexed
            
            # Create mutant sequence
            mutant_seq = list(reference_sequence)
            mutant_seq[mutation_pos] = pred['to']
            mutant_seq = ''.join(mutant_seq)
            
            # Predict mutant structure
            print(f"Validating {pred['mutation']}...")
            try:
                mut_structure = self.neurosnap.predict_structure(
                    mutant_seq,
                    job_name=f"D1_{pred['mutation']}"
                )
                mut_plddt = mut_structure['plddt_scores']
                
                # Calculate confidence change at mutation site (±5 residues)
                window_start = max(0, mutation_pos - 5)
                window_end = min(len(ref_plddt), mutation_pos + 6)
                
                ref_confidence = np.mean(ref_plddt[window_start:window_end])
                mut_confidence = np.mean(mut_plddt[window_start:window_end])
                confidence_change = mut_confidence - ref_confidence
                
                structure_confidence_changes.append(confidence_change)
                predicted_ddgs.append(pred['delta_stability'])
                
                validated_predictions.append({
                    **pred,
                    'structure_confidence_change': confidence_change,
                    'alphafold_confidence': mut_confidence,
                    'structure_validated': True
                })
                
            except Exception as e:
                print(f"Structure prediction failed for {pred['mutation']}: {e}")
                validated_predictions.append({
                    **pred,
                    'structure_validated': False
                })
        
        # Calculate correlation between predicted ΔΔG and structure confidence changes
        if len(structure_confidence_changes) > 3:
            # Negative correlation expected: more negative ΔΔG → higher confidence
            r_value, p_value = stats.pearsonr(predicted_ddgs, structure_confidence_changes)
            r_squared = r_value ** 2
            
            # Also calculate Spearman for non-linear relationships
            spearman_r, spearman_p = stats.spearmanr(predicted_ddgs, structure_confidence_changes)
            
            return {
                'pearson_r': r_value,
                'pearson_r2': r_squared,
                'pearson_p_value': p_value,
                'spearman_r': spearman_r,
                'spearman_p_value': spearman_p,
                'n_validated': len(structure_confidence_changes),
                'validated_predictions': validated_predictions,
                'interpretation': self._interpret_correlation(r_squared, p_value)
            }
        else:
            return {
                'error': 'Insufficient data for correlation',
                'validated_predictions': validated_predictions
            }
    
    def _interpret_correlation(self, r2: float, p_value: float) -> str:
        """Interpret correlation strength and significance"""
        if p_value > 0.05:
            return f"No significant correlation (p={p_value:.3f})"
        
        strength = "weak" if r2 < 0.3 else "moderate" if r2 < 0.7 else "strong"
        return f"{strength.capitalize()} significant correlation (R²={r2:.3f}, p={p_value:.3f})"
    
    def benchmark_with_neurofold(self,
                                 reference_sequence: str,
                                 local_predictions: List[Dict]) -> Dict:
        """
        Benchmark local ML predictions against NeuroFold's predictions
        """
        print("Running NeuroFold benchmark...")
        
        # Get NeuroFold's optimization suggestions
        neurofold_variants = self.neurosnap.optimize_thermostability(
            reference_sequence,
            num_designs=10
        )
        
        # Compare predicted ΔTm between local and NeuroFold
        local_delta_tms = [pred['estimated_delta_tm'] for pred in local_predictions[:10]]
        neurofold_delta_tms = [var['predicted_tm_increase'] for var in neurofold_variants]
        
        # Calculate agreement statistics
        if len(local_delta_tms) == len(neurofold_delta_tms):
            correlation_r, correlation_p = stats.pearsonr(local_delta_tms, neurofold_delta_tms)
            
            return {
                'correlation_r': correlation_r,
                'correlation_r2': correlation_r ** 2,
                'p_value': correlation_p,
                'mean_absolute_error': np.mean(np.abs(np.array(local_delta_tms) - np.array(neurofold_delta_tms))),
                'neurofold_predictions': neurofold_variants,
                'agreement_interpretation': self._interpret_correlation(correlation_r ** 2, correlation_p)
            }
        
        return {'error': 'Mismatched prediction counts'}
    
    def cross_species_structure_validation(self,
                                          sequences: Dict[str, str]) -> Dict:
        """
        Validate conservation analysis using AlphaFold2 structures across species
        """
        print("Running cross-species structure validation...")
        
        structure_similarities = []
        species_pairs = []
        
        # Predict structures for each species
        structures = {}
        for species, sequence in list(sequences.items())[:5]:  # Limit to 5 species
            try:
                print(f"Predicting structure for {species}...")
                structure = self.neurosnap.predict_structure(sequence, job_name=f"D1_{species}")
                structures[species] = structure
            except Exception as e:
                print(f"Failed for {species}: {e}")
        
        # Calculate structural similarity (using pLDDT correlation as proxy)
        species_list = list(structures.keys())
        for i in range(len(species_list)):
            for j in range(i + 1, len(species_list)):
                sp1, sp2 = species_list[i], species_list[j]
                
                plddt1 = structures[sp1]['plddt_scores']
                plddt2 = structures[sp2]['plddt_scores']
                
                # Align to same length
                min_len = min(len(plddt1), len(plddt2))
                plddt1 = plddt1[:min_len]
                plddt2 = plddt2[:min_len]
                
                # Calculate correlation
                r, p = stats.pearsonr(plddt1, plddt2)
                
                structure_similarities.append(r)
                species_pairs.append(f"{sp1}-{sp2}")
        
        return {
            'species_pairs': species_pairs,
            'structure_similarities': structure_similarities,
            'mean_similarity': np.mean(structure_similarities),
            'std_similarity': np.std(structure_similarities),
            'structures': structures
        }


# Example usage integration with existing pipeline
class HybridPredictionPipeline:
    """
    Hybrid pipeline combining local ML models with NeuroSnap validation
    """
    
    def __init__(self, neurosnap_api_key: Optional[str] = None):
        from ml_models import MutationEffectPredictor, VariantScorer
        
        self.local_predictor = MutationEffectPredictor()
        self.local_scorer = VariantScorer()
        
        if neurosnap_api_key:
            self.neurosnap = NeuroSnapAPI(neurosnap_api_key)
            self.analyzer = EnhancedStatisticalAnalyzer(self.neurosnap)
            self.use_neurosnap = True
        else:
            self.use_neurosnap = False
            print("Warning: No NeuroSnap API key provided. Running in local-only mode.")
    
    def analyze_with_validation(self,
                               reference_sequence: str,
                               species_sequences: Dict[str, str]) -> Dict:
        """
        Complete analysis with NeuroSnap validation
        """
        from analysis_engine import run_comprehensive_analysis
        
        # Run local analysis
        print("Running local comparative analysis...")
        local_results = run_comprehensive_analysis(
            list(species_sequences.keys()),
            reference_sequence
        )
        
        if not self.use_neurosnap:
            return {
                **local_results,
                'neurosnap_validation': None,
                'mode': 'local_only'
            }
        
        # Enhance with NeuroSnap validation
        print("\nEnhancing with NeuroSnap validation...")
        
        # Validate top predictions with structure
        structure_validation = self.analyzer.validate_predictions_with_structure(
            local_results['recommendations'],
            reference_sequence
        )
        
        # Benchmark against NeuroFold
        neurofold_benchmark = self.analyzer.benchmark_with_neurofold(
            reference_sequence,
            local_results['recommendations']
        )
        
        # Cross-species structure validation
        structure_comparison = self.analyzer.cross_species_structure_validation(
            species_sequences
        )
        
        return {
            **local_results,
            'neurosnap_validation': {
                'structure_validation': structure_validation,
                'neurofold_benchmark': neurofold_benchmark,
                'cross_species_structures': structure_comparison
            },
            'mode': 'hybrid_validated',
            'enhanced_statistics': {
                'local_vs_alphafold_r2': structure_validation.get('pearson_r2', 0),
                'local_vs_alphafold_p': structure_validation.get('pearson_p_value', 1),
                'local_vs_neurofold_r2': neurofold_benchmark.get('correlation_r2', 0),
                'local_vs_neurofold_p': neurofold_benchmark.get('p_value', 1),
                'mean_structure_similarity': structure_comparison.get('mean_similarity', 0)
            }
        }