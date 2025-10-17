"""
neurosnap_client.py
Proper NeuroSnap API client for D1 protein engineering
Replaces the mock implementation with real API calls
"""

import aiohttp
import asyncio
import os
from typing import List, Dict, Optional
import json
from pathlib import Path


class NeurosnapClient:
    """
    Production-ready NeuroSnap API client
    
    Available tools:
    - TemStaPro: Thermostability prediction (40-65°C)
    - ThermoMPNN: Structure-based mutation stability
    - AlphaFold2: High-accuracy structure prediction
    - RFdiffusion: AI-powered protein design
    - ESM-IF1: Inverse folding for sequence design
    """
    
    def __init__(self, api_key: str):
        self.api_key = api_key
        self.base_url = "https://api.neurosnap.ai/v1"
        self.headers = {
            "Authorization": f"Bearer {api_key}",
            "Content-Type": "application/json"
        }
        self.session = None
    
    async def __aenter__(self):
        """Async context manager entry"""
        self.session = aiohttp.ClientSession()
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit"""
        if self.session:
            await self.session.close()
    
    # ==================== TEMSTAPRO: Thermostability Prediction ====================
    
    async def predict_thermostability_temstapro(self,
                                                sequence: str,
                                                temperature: int = 50,
                                                job_name: str = "D1_thermostability") -> Dict:
        """
        Predict if protein is stable at given temperature using TemStaPro
        
        Args:
            sequence: Protein sequence (single letter amino acids)
            temperature: Target temperature in Celsius (40-65°C)
            job_name: Name for tracking the job
        
        Returns:
            {
                'job_id': str,
                'status': 'completed' | 'running' | 'failed',
                'prediction': {
                    'is_stable': bool,
                    'confidence': float (0-1),
                    'predicted_tm': float  # melting temperature
                }
            }
        """
        if not self.session:
            self.session = aiohttp.ClientSession()
        
        payload = {
            "sequence": sequence,
            "temperature_threshold": temperature,
            "job_name": job_name
        }
        
        try:
            # Submit job
            async with self.session.post(
                f"{self.base_url}/temstapro/submit",
                headers=self.headers,
                json=payload
            ) as response:
                response.raise_for_status()
                result = await response.json()
                job_id = result.get('job_id')
            
            # Poll for results
            return await self._wait_for_job(job_id, tool_name="temstapro")
        
        except aiohttp.ClientError as e:
            print(f"❌ TemStaPro API error: {e}")
            return self._fallback_thermostability(sequence, temperature)
    
    async def batch_predict_temstapro(self,
                                     sequences: List[str],
                                     temperature: int = 50) -> List[Dict]:
        """Batch prediction for multiple sequences"""
        tasks = []
        for i, seq in enumerate(sequences):
            task = self.predict_thermostability_temstapro(
                seq, 
                temperature, 
                job_name=f"D1_variant_{i}"
            )
            tasks.append(task)
        
        return await asyncio.gather(*tasks)
    
    # ==================== THERMOMPNN: Mutation Stability ====================
    
    async def predict_mutation_effects_thermompnn(self,
                                                 structure_pdb: str,
                                                 mutations: List[Dict],
                                                 job_name: str = "D1_mutations") -> List[Dict]:
        """
        Predict stability changes for mutations using ThermoMPNN
        
        Args:
            structure_pdb: PDB file content or path to PDB
            mutations: List of mutations, e.g. [{'position': 123, 'from': 'A', 'to': 'V'}]
            job_name: Job identifier
        
        Returns:
            List of mutation effects with ΔΔG predictions
            [{
                'mutation': 'A123V',
                'predicted_ddg': float,  # kcal/mol (negative = stabilizing)
                'confidence': float (0-1),
                'recommendation': 'stabilizing' | 'neutral' | 'destabilizing'
            }]
        """
        if not self.session:
            self.session = aiohttp.ClientSession()
        
        # Format mutations for ThermoMPNN
        formatted_mutations = [
            f"{m['from']}{m['position']}{m['to']}" 
            for m in mutations
        ]
        
        payload = {
            "structure": structure_pdb,
            "mutations": formatted_mutations,
            "job_name": job_name
        }
        
        try:
            async with self.session.post(
                f"{self.base_url}/thermompnn/submit",
                headers=self.headers,
                json=payload
            ) as response:
                response.raise_for_status()
                result = await response.json()
                job_id = result.get('job_id')
            
            return await self._wait_for_job(job_id, tool_name="thermompnn")
        
        except aiohttp.ClientError as e:
            print(f"❌ ThermoMPNN API error: {e}")
            return self._fallback_mutation_prediction(mutations)
    
    # ==================== ALPHAFOLD2: Structure Prediction ====================
    
    async def predict_structure_alphafold(self,
                                         sequence: str,
                                         job_name: str = "D1_structure",
                                         use_templates: bool = False) -> Dict:
        """
        Predict 3D structure using AlphaFold2
        
        Args:
            sequence: Protein sequence
            job_name: Job identifier
            use_templates: Whether to use template structures (slower but more accurate)
        
        Returns:
            {
                'job_id': str,
                'pdb_content': str,  # PDB file content
                'confidence': float,  # pLDDT score (0-100)
                'pdb_url': str  # URL to download PDB
            }
        """
        if not self.session:
            self.session = aiohttp.ClientSession()
        
        payload = {
            "sequence": sequence,
            "job_name": job_name,
            "use_templates": use_templates,
            "num_models": 1  # Generate 1 model (faster)
        }
        
        try:
            async with self.session.post(
                f"{self.base_url}/alphafold2/submit",
                headers=self.headers,
                json=payload
            ) as response:
                response.raise_for_status()
                result = await response.json()
                job_id = result.get('job_id')
            
            print(f"📊 AlphaFold2 job submitted: {job_id}")
            print(f"⏳ This may take 5-15 minutes depending on sequence length...")
            
            return await self._wait_for_job(job_id, tool_name="alphafold2")
        
        except aiohttp.ClientError as e:
            print(f"❌ AlphaFold2 API error: {e}")
            return {"error": str(e), "pdb_content": None}
    
    # ==================== RFDIFFUSION: Protein Design ====================
    
    async def generate_variants_rfdiffusion(self,
                                           template_pdb: str,
                                           design_regions: List[tuple],
                                           num_designs: int = 50,
                                           job_name: str = "D1_design") -> List[Dict]:
        """
        Generate novel protein variants using RFdiffusion
        
        Args:
            template_pdb: Starting structure (PDB content)
            design_regions: List of (start, end) tuples for regions to redesign
                           e.g. [(50, 100), (150, 200)]
            num_designs: Number of variants to generate
            job_name: Job identifier
        
        Returns:
            List of designed sequences with confidence scores
        """
        if not self.session:
            self.session = aiohttp.ClientSession()
        
        payload = {
            "template_structure": template_pdb,
            "design_regions": design_regions,
            "num_designs": num_designs,
            "job_name": job_name
        }
        
        try:
            async with self.session.post(
                f"{self.base_url}/rfdiffusion/submit",
                headers=self.headers,
                json=payload
            ) as response:
                response.raise_for_status()
                result = await response.json()
                job_id = result.get('job_id')
            
            print(f"🧬 RFdiffusion job submitted: {job_id}")
            print(f"⏳ Generating {num_designs} designs, this may take 10-30 minutes...")
            
            return await self._wait_for_job(job_id, tool_name="rfdiffusion")
        
        except aiohttp.ClientError as e:
            print(f"❌ RFdiffusion API error: {e}")
            return []
    
    # ==================== HELPER METHODS ====================
    
    async def _wait_for_job(self, 
                           job_id: str, 
                           tool_name: str,
                           max_wait_seconds: int = 1800,
                           poll_interval: int = 10) -> Dict:
        """
        Poll job status until completion
        
        Args:
            job_id: Job identifier from submission
            tool_name: Name of the tool (for endpoint)
            max_wait_seconds: Maximum time to wait (default 30 min)
            poll_interval: Seconds between status checks
        """
        if not self.session:
            self.session = aiohttp.ClientSession()
        
        elapsed = 0
        while elapsed < max_wait_seconds:
            try:
                async with self.session.get(
                    f"{self.base_url}/{tool_name}/status/{job_id}",
                    headers=self.headers
                ) as response:
                    response.raise_for_status()
                    status = await response.json()
                
                if status['status'] == 'completed':
                    print(f"✅ {tool_name} job {job_id} completed!")
                    return status['results']
                elif status['status'] == 'failed':
                    print(f"❌ {tool_name} job {job_id} failed: {status.get('error')}")
                    return {"error": status.get('error')}
                else:
                    # Still running
                    progress = status.get('progress', 0)
                    print(f"⏳ {tool_name} progress: {progress}%")
                    await asyncio.sleep(poll_interval)
                    elapsed += poll_interval
            
            except aiohttp.ClientError as e:
                print(f"⚠️  Error checking job status: {e}")
                await asyncio.sleep(poll_interval)
                elapsed += poll_interval
        
        print(f"⏰ Job {job_id} timed out after {max_wait_seconds}s")
        return {"error": "Job timed out"}
    
    def _fallback_thermostability(self, sequence: str, temperature: int) -> Dict:
        """Fallback heuristic when API fails"""
        # Simple rule: count stabilizing residues
        hydrophobic = sum(aa in 'ILVAFW' for aa in sequence) / len(sequence)
        charged = sum(aa in 'KRED' for aa in sequence) / len(sequence)
        
        stability_score = hydrophobic * 0.6 + charged * 0.4
        is_stable = stability_score > 0.5
        
        return {
            'prediction': {
                'is_stable': is_stable,
                'confidence': 0.3,  # Low confidence for fallback
                'predicted_tm': 40 + (stability_score * 30)
            },
            'note': 'Using fallback heuristic - API unavailable'
        }
    
    def _fallback_mutation_prediction(self, mutations: List[Dict]) -> List[Dict]:
        """Fallback when ThermoMPNN fails"""
        results = []
        for mut in mutations:
            # Simple heuristic: hydrophobic substitutions stabilizing
            from_aa = mut['from']
            to_aa = mut['to']
            
            hydrophobic = set('ILVAFW')
            if to_aa in hydrophobic and from_aa not in hydrophobic:
                ddg = -1.5  # Stabilizing
            elif from_aa in hydrophobic and to_aa not in hydrophobic:
                ddg = 1.5  # Destabilizing
            else:
                ddg = 0.0  # Neutral
            
            results.append({
                'mutation': f"{from_aa}{mut['position']}{to_aa}",
                'predicted_ddg': ddg,
                'confidence': 0.3,
                'recommendation': 'stabilizing' if ddg < 0 else 'destabilizing' if ddg > 0 else 'neutral',
                'note': 'Fallback prediction - API unavailable'
            })
        
        return results
    
    # ==================== UTILITY METHODS ====================
    
    async def get_account_info(self) -> Dict:
        """Get account information and credits"""
        if not self.session:
            self.session = aiohttp.ClientSession()
        
        try:
            async with self.session.get(
                f"{self.base_url}/account",
                headers=self.headers
            ) as response:
                response.raise_for_status()
                return await response.json()
        except aiohttp.ClientError as e:
            print(f"Error fetching account info: {e}")
            return {"error": str(e)}
    
    async def list_jobs(self, limit: int = 10) -> List[Dict]:
        """List recent jobs"""
        if not self.session:
            self.session = aiohttp.ClientSession()
        
        try:
            async with self.session.get(
                f"{self.base_url}/jobs?limit={limit}",
                headers=self.headers
            ) as response:
                response.raise_for_status()
                return await response.json()
        except aiohttp.ClientError as e:
            print(f"Error listing jobs: {e}")
            return []


# ==================== CONVENIENCE FUNCTIONS ====================

async def quick_thermostability_check(sequence: str, 
                                     target_temp: int = 50,
                                     api_key: Optional[str] = None) -> Dict:
    """
    Quick function to check thermostability of a single sequence
    
    Usage:
        result = await quick_thermostability_check("MTVK...", target_temp=55)
    """
    if api_key is None:
        api_key = os.getenv('NEUROSNAP_API_KEY')
    
    async with NeurosnapClient(api_key) as client:
        return await client.predict_thermostability_temstapro(sequence, target_temp)


async def quick_structure_prediction(sequence: str,
                                    api_key: Optional[str] = None) -> Dict:
    """
    Quick function to predict structure
    
    Usage:
        result = await quick_structure_prediction("MTVK...")
    """
    if api_key is None:
        api_key = os.getenv('NEUROSNAP_API_KEY')
    
    async with NeurosnapClient(api_key) as client:
        return await client.predict_structure_alphafold(sequence)


# ==================== TESTING ====================

async def test_neurosnap_connection():
    """Test connection to NeuroSnap API"""
    api_key = os.getenv('NEUROSNAP_API_KEY')
    
    if not api_key:
        print("❌ No API key found in environment")
        print("   Set NEUROSNAP_API_KEY in your .env file")
        return False
    
    print("🔍 Testing NeuroSnap API connection...")
    
    async with NeurosnapClient(api_key) as client:
        # Test account access
        account = await client.get_account_info()
        
        if 'error' in account:
            print(f"❌ Connection failed: {account['error']}")
            return False
        
        print(f"✅ Connected successfully!")
        print(f"   Account: {account.get('email', 'N/A')}")
        print(f"   Credits: {account.get('credits', 'N/A')}")
        
        return True


if __name__ == "__main__":
    # Test the client
    asyncio.run(test_neurosnap_connection())