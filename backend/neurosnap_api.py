"""
neurosnap_api.py
CRITICAL FILE - Required by master_pipeline_v1.py
NeuroSnap API integration for AlphaFold2 and NeuroFold validation
"""

import requests
import time
from typing import Dict, List, Optional
from requests_toolbelt.multipart.encoder import MultipartEncoder
import numpy as np
from scipy import stats


class NeuroSnapAPI:
    """NeuroSnap API client for protein structure prediction and validation"""
    
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
            data=multipart_data,
            timeout=30
        )
        
        if response.status_code == 200:
            return response.json()
        else:
            raise Exception(f"Job submission failed: {response.text}")
    
    def _check_status(self, job_id: str) -> str:
        """Check job status"""
        response = requests.get(
            f"{self.base_url}/job/status/{job_id}",
            headers=self.headers,
            timeout=10
        )
        return response.json()
    
    def _get_results(self, job_id: str) -> Dict:
        """Get job results"""
        response = requests.get(
            f"{self.base_url}/job/output/{job_id}",
            headers=self.headers,
            timeout=30
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
            
            time.sleep(10)
        
        raise TimeoutError(f"Job {job_id} timed out")
    
    def predict_structure(self, 
                         sequence: str, 
                         job_name: str = "D1_structure",
                         use_templates: bool = False) -> Dict:
        """Predict protein structure using AlphaFold2"""
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
            'pdb_structure': results.get('ranked_0.pdb', ''),
            'plddt_scores': self._parse_plddt(results),
            'confidence': np.mean(self._parse_plddt(results)) if self._parse_plddt(results).size > 0 else 0.0
        }
    
    def _parse_plddt(self, results: Dict) -> np.ndarray:
        """Extract pLDDT scores from AlphaFold2 results"""
        pdb_content = results.get('ranked_0.pdb', '')
        plddt_scores = []
        
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                try:
                    plddt = float(line[60:66].strip())
                    plddt_scores.append(plddt)
                except (ValueError, IndexError):
                    continue
        
        return np.array(plddt_scores)
    
    def optimize_thermostability(self,
                                 sequence: str,
                                 num_designs: int = 10,
                                 job_name: str = "D1_thermostability") -> List[Dict]:
        """Optimize protein for thermostability using NeuroFold"""
        fields = {
            "Sequence": sequence,
            "Number of Designs": str(num_designs),
            "Optimization Target": "Thermostability",
            "Note": job_name
        }
        
        job_id = self._submit_job("NeuroFold", fields)
        results = self._wait_for_completion(job_id, timeout=7200)
        
        variants = []
        for i in range(num_designs):
            variant_key = f"design_{i}"
            if variant_key in results:
                variants.append({
                    'sequence': results[variant_key]['sequence'],
                    'predicted_tm_increase': results[variant_key].get('delta_tm', 0),
                    'confidence': results[variant_key].get('confidence', 0)
                })
        
        return variants


class EnhancedStatisticalAnalyzer:
    """Statistical analysis combining local predictions with NeuroSnap validation"""
    
    def __init__(self, neurosnap_api: NeuroSnapAPI):
        self.neurosnap = neurosnap_api
    
    def validate_predictions_with_structure(self,
                                           predictions: List[Dict],
                                           reference_sequence: str) -> Dict:
        """Validate mutation predictions using AlphaFold2 structures"""
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
        
        # Validate top 10 mutations
        for pred in predictions[:10]:
            mutation_pos = pred['position'] - 1
            
            # Create mutant sequence
            mutant_seq = list(reference_sequence)
            mutant_seq[mutation_pos] = pred['to']
            mutant_seq = ''.join(mutant_seq)
            
            print(f"Validating {pred['mutation']}...")
            try:
                mut_structure = self.neurosnap.predict_structure(
                    mutant_seq,
                    job_name=f"D1_{pred['mutation']}"
                )
                mut_plddt = mut_structure['plddt_scores']
                
                # Calculate confidence change
                window_start = max(0, mutation_pos - 5)
                window_end = min(len(ref_plddt), mutation_pos + 6)
                
                ref_confidence = np.mean(ref_plddt[window_start:window_end])
                mut_confidence = np.mean(mut_plddt[window_start:window_end])
                confidence_change = mut_confidence - ref_confidence
                
                structure_confidence_changes.append(confidence_change)
                predicted_ddgs.append(pred.get('delta_stability', 0))
                
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
        
        # Calculate correlation
        if len(structure_confidence_changes) > 3:
            r_value, p_value = stats.pearsonr(predicted_ddgs, structure_confidence_changes)
            r_squared = r_value ** 2
            spearman_r, spearman_p = stats.spearmanr(predicted_ddgs, structure_confidence_changes)
            
            return {
                'pearson_r': r_value,
                'pearson_r2': r_squared,
                'pearson_p_value': p_value,
                'spearman_r': spearman_r,
                'spearman_p_value': spearman_p,
                'n_validated': len(structure_confidence_changes),
                'validated_predictions': validated_predictions
            }
        
        return {
            'error': 'Insufficient data',
            'validated_predictions': validated_predictions
        }


class HybridPredictionPipeline:
    """Hybrid pipeline combining local ML models with NeuroSnap validation"""
    
    def __init__(self, neurosnap_api_key: Optional[str] = None):
        if neurosnap_api_key:
            self.neurosnap = NeuroSnapAPI(neurosnap_api_key)
            self.analyzer = EnhancedStatisticalAnalyzer(self.neurosnap)
            self.use_neurosnap = True
            print("✓ NeuroSnap API initialized")
        else:
            self.use_neurosnap = False
            print("⚠ No NeuroSnap API key - running in local-only mode")
    
    def analyze_with_validation(self,
                               reference_sequence: str,
                               species_sequences: Dict[str, str]) -> Dict:
        """Complete analysis with NeuroSnap validation"""
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
        
        structure_validation = self.analyzer.validate_predictions_with_structure(
            local_results['recommendations'],
            reference_sequence
        )
        
        return {
            **local_results,
            'neurosnap_validation': {
                'structure_validation': structure_validation
            },
            'mode': 'hybrid_validated',
            'enhanced_statistics': {
                'local_vs_alphafold_r2': structure_validation.get('pearson_r2', 0),
                'local_vs_alphafold_p': structure_validation.get('pearson_p_value', 1)
            }
        }