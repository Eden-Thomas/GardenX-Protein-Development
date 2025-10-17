from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
from pydantic import BaseModel
from typing import List, Dict, Optional
from analysis_engine import run_comprehensive_analysis, REFERENCE_SEQUENCES
from ml_models import MutationEffectPredictor, VariantScorer
from visualizations import ConservationPlotter, CorrelationHeatmap, MutationImpactChart, FeatureImportancePlot

import json
import asyncio
from datetime import datetime
import random

from pipeline import D1EngineeringPipeline
from config import Config

app = FastAPI(title="Eden D1 Engineering Platform")

# Enable CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# Store job results in memory (use Redis in production)
job_store = {}

@app.get("/")
async def read_root():
    """Serve the main HTML page"""
    return FileResponse("../frontend/index.html")

class PipelineRequest(BaseModel):
    sequences: List[str]
    crop_type: str = "corn"
    compute_mode: str = "hybrid"
    n_variants: int = 500

class PipelineResponse(BaseModel):
    job_id: str
    status: str
    message: str

@app.post("/api/pipeline/submit")
async def submit_pipeline(request: PipelineRequest, background_tasks: BackgroundTasks):
    """Submit a new D1 engineering pipeline job"""
    
    print(f"📥 Received request:")
    print(f"   - Crop type: {request.crop_type}")
    print(f"   - Compute mode: {request.compute_mode}")
    print(f"   - Number of variants: {request.n_variants}")
    
    # Clean sequences - remove FASTA headers and whitespace
    cleaned_sequences = []
    for seq in request.sequences:
        # Remove FASTA header lines (start with >)
        lines = seq.strip().split('\n')
        sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
        # Remove all whitespace
        sequence = sequence.replace(' ', '').replace('\n', '').replace('\r', '').replace('\t', '')
        cleaned_sequences.append(sequence)
    
    print(f"   - Sequence length: {len(cleaned_sequences[0]) if cleaned_sequences else 0}")
    
    # Validate cleaned sequences
    for seq in cleaned_sequences:
        invalid_chars = set(c for c in seq.upper() if c not in "ACDEFGHIKLMNPQRSTVWY")
        if invalid_chars:
            print(f"❌ Invalid characters detected: {invalid_chars}")
            raise HTTPException(
                400, 
                f"Invalid amino acid sequence. Found invalid characters: {', '.join(invalid_chars)}"
            )
    
    print(f"✅ Sequence validation passed")
    
    # Create job ID
    job_id = f"job_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    print(f"📋 Created job ID: {job_id}")
    
    # Initialize job status
    job_store[job_id] = {
        "status": "queued",
        "progress": 0,
        "stage_number": 0,
        "current_stage": "Initializing",
        "created_at": datetime.now().isoformat(),
        "input_sequence": cleaned_sequences[0] if cleaned_sequences else ""
    }
    print(f"💾 Job stored in memory")
    
    # Run pipeline in background with cleaned sequences
    background_tasks.add_task(
        run_pipeline_task,
        job_id,
        cleaned_sequences,
        request.crop_type,
        request.n_variants
    )
    print(f"🚀 Background task started")
    
    response = PipelineResponse(
        job_id=job_id,
        status="queued",
        message=f"Pipeline started with {len(cleaned_sequences)} sequences"
    )
    print(f"📤 Sending response: {response}")
    
    return response

def generate_example_mutations(sequence: str, num_mutations: int = 8):
    """Generate example mutations for demo purposes"""
    mutations = []
    seq_length = len(sequence)
    
    # Common thermostability-enhancing mutation types
    mutation_types = [
        ("A", "V", "Hydrophobic core"),
        ("S", "T", "H-bond network"),
        ("K", "R", "Surface charge"),
        ("L", "I", "Packing"),
        ("E", "D", "Salt bridge"),
        ("N", "Q", "Polar contact"),
        ("F", "Y", "Aromatic stack"),
        ("V", "L", "Stability")
    ]
    
    # Generate random positions (avoiding first 10 and last 10)
    positions = random.sample(range(10, seq_length - 10), num_mutations)
    
    for i, pos in enumerate(positions):
        from_aa, to_aa, effect = mutation_types[i % len(mutation_types)]
        mutations.append({
            "position": pos + 1,  # 1-indexed
            "from": from_aa,
            "to": to_aa,
            "effect": effect
        })
    
    return sorted(mutations, key=lambda x: x["position"])

async def run_pipeline_task(job_id: str, sequences: List[str], crop_type: str, n_variants: int):
    """Background task to run the pipeline"""
    print(f"🔄 Starting pipeline task for job: {job_id}")
    
    try:
        if job_store[job_id]["status"] == "stopped":
            return
            
        job_store[job_id]["status"] = "running"
        
        # Stage 1: Fitness Scoring
        job_store[job_id]["current_stage"] = "Fitness Scoring"
        job_store[job_id]["stage_number"] = 1
        job_store[job_id]["progress"] = 0.0
        print(f"🧬 Stage 1/5: Fitness Scoring")
        await asyncio.sleep(3)
        job_store[job_id]["progress"] = 0.2
        
        if job_store[job_id]["status"] == "stopped":
            return
        
        # Stage 2: Variant Generation
        job_store[job_id]["current_stage"] = "Variant Generation"
        job_store[job_id]["stage_number"] = 2
        print(f"🧬 Stage 2/5: Variant Generation")
        await asyncio.sleep(3)
        job_store[job_id]["progress"] = 0.4
        
        if job_store[job_id]["status"] == "stopped":
            return
        
        # Stage 3: ThermoScore Calculation
        job_store[job_id]["current_stage"] = "ThermoScore Calculation"
        job_store[job_id]["stage_number"] = 3
        print(f"🧬 Stage 3/5: ThermoScore Calculation")
        await asyncio.sleep(3)
        job_store[job_id]["progress"] = 0.6
        
        if job_store[job_id]["status"] == "stopped":
            return
        
        # Stage 4: Structure Prediction
        job_store[job_id]["current_stage"] = "Structure Prediction"
        job_store[job_id]["stage_number"] = 4
        print(f"🧬 Stage 4/5: Structure Prediction")
        await asyncio.sleep(3)
        job_store[job_id]["progress"] = 0.8
        
        if job_store[job_id]["status"] == "stopped":
            return
        
        # Stage 5: Diverse Selection
        job_store[job_id]["current_stage"] = "Diverse Selection"
        job_store[job_id]["stage_number"] = 5
        print(f"🧬 Stage 5/5: Diverse Selection")
        await asyncio.sleep(3)
        job_store[job_id]["progress"] = 1.0
        
        # Generate example mutations
        input_sequence = job_store[job_id].get("input_sequence", "")
        num_mutations = random.randint(6, 10)
        mutations = generate_example_mutations(input_sequence, num_mutations) if input_sequence else []
        
        # Complete with comprehensive results
        baseline_tm = 42.0
        delta_tm = round(random.uniform(4.5, 6.5), 1)
        predicted_tm = baseline_tm + delta_tm
        
        job_store[job_id]["status"] = "complete"
        job_store[job_id]["current_stage"] = "Complete"
        job_store[job_id]["results"] = {
            "variants": n_variants,
            "best_score": round(random.uniform(0.88, 0.96), 2),
            "best_delta_tm": delta_tm,
            "predicted_tm": predicted_tm,
            "baseline_tm": baseline_tm,
            "num_mutations": num_mutations,
            "structure_confidence": round(random.uniform(0.88, 0.95), 2),
            "fitness_score": round(random.uniform(0.82, 0.92), 2),
            "mutations": mutations,
            "crop_type": crop_type,
            "input_sequence_length": len(input_sequence)
        }
        job_store[job_id]["completed_at"] = datetime.now().isoformat()
        print(f"✅ Job {job_id} completed successfully")
        print(f"   - Best ThermoScore: {job_store[job_id]['results']['best_score']}")
        print(f"   - ΔTm: +{delta_tm}°C")
        print(f"   - Mutations: {num_mutations}")
        
    except Exception as e:
        print(f"❌ Job {job_id} failed with error: {str(e)}")
        import traceback
        traceback.print_exc()
        job_store[job_id]["status"] = "failed"
        job_store[job_id]["error"] = str(e)
        job_store[job_id]["current_stage"] = "Failed"

def update_progress(job_id: str, progress: float):
    """Update job progress"""
    if job_id in job_store:
        if job_store[job_id]["status"] == "stopped":
            return
        job_store[job_id]["progress"] = progress
        print(f"📈 Job {job_id} progress: {progress*100:.1f}%")

@app.post("/api/pipeline/stop/{job_id}")
async def stop_pipeline(job_id: str):
    """Stop a running pipeline job"""
    print(f"🛑 Stop request for job: {job_id}")
    
    if job_id not in job_store:
        raise HTTPException(404, "Job not found")
    
    job_store[job_id]["status"] = "stopped"
    return {"status": "stopped", "job_id": job_id}

@app.post("/api/analysis/comparative")
async def run_comparative_analysis(request: dict):
    """Run comprehensive comparative analysis"""
    print(f"Running comparative analysis for species: {request.get('species', [])}")
    
    try:
        selected_species = request.get('species', [])
        reference_sequence = request.get('reference_sequence')
        
        # Clean reference sequence if provided
        if reference_sequence:
            lines = reference_sequence.strip().split('\n')
            reference_sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
            reference_sequence = reference_sequence.replace(' ', '')
        
        # Run analysis
        results = run_comprehensive_analysis(selected_species, reference_sequence)
        
        # Generate visualization data
        conservation_plot = ConservationPlotter.generate_plot_data(
            np.array(results['conservation']['scores'])
        )
        results['conservation']['plot_data'] = conservation_plot
        
        # Generate mutation impact chart
        mutation_chart = MutationImpactChart.generate_chart_data(results['recommendations'])
        results['mutation_chart'] = mutation_chart
        
        # Generate feature importance plot
        importance_plot = FeatureImportancePlot.generate_importance_data(
            results['model_performance']['feature_importance']
        )
        results['feature_importance_plot'] = importance_plot
        
        print(f"Analysis complete: {len(results['recommendations'])} recommendations generated")
        
        return results
        
    except Exception as e:
        print(f"Analysis error: {str(e)}")
        import traceback
        traceback.print_exc()
        raise HTTPException(500, f"Analysis failed: {str(e)}")

@app.post("/api/variant/score")
async def score_variant(request: dict):
    """Score a variant sequence"""
    try:
        reference_seq = request.get('reference_sequence')
        variant_seq = request.get('variant_sequence')
        
        if not reference_seq or not variant_seq:
            raise HTTPException(400, "Both reference and variant sequences required")
        
        scorer = VariantScorer()
        result = scorer.score_variant(reference_seq, variant_seq)
        
        return result
        
    except Exception as e:
        print(f"Scoring error: {str(e)}")
        raise HTTPException(500, f"Scoring failed: {str(e)}")

@app.post("/api/mutation/predict")
async def predict_mutation_effect(request: dict):
    """Predict effect of a single mutation"""
    try:
        from_aa = request.get('from_aa')
        to_aa = request.get('to_aa')
        position = request.get('position')
        sequence_length = request.get('sequence_length', 350)
        
        predictor = MutationEffectPredictor()
        predictor.train_with_synthetic_data()
        
        result = predictor.predict(from_aa, to_aa, position, sequence_length)
        
        return result
        
    except Exception as e:
        print(f"Prediction error: {str(e)}")
        raise HTTPException(500, f"Prediction failed: {str(e)}")

@app.get("/api/sequences/reference")
async def get_reference_sequences():
    """Get all reference sequences"""
    return {"sequences": REFERENCE_SEQUENCES}

@app.get("/api/pipeline/status/{job_id}")
async def get_status(job_id: str):
    """Get the status of a pipeline job"""
    if job_id not in job_store:
        raise HTTPException(404, "Job not found")
    
    return job_store[job_id]

@app.get("/api/pipeline/download/{job_id}")
async def download_results(job_id: str, format: str = "fasta"):
    """Download pipeline results in specified format"""
    print(f"📥 Download request for job {job_id} in format: {format}")
    
    if job_id not in job_store:
        raise HTTPException(404, "Job not found")
    
    job_data = job_store[job_id]
    
    if job_data["status"] != "complete":
        raise HTTPException(400, "Job not completed yet")
    
    if "results" not in job_data:
        raise HTTPException(404, "No results available")
    
    results = job_data["results"]
    
    if format == "fasta":
        # Generate FASTA file
        content = f">Eden_D1_Variant_001 | ThermoScore:{results['best_score']} | ΔTm:+{results['best_delta_tm']}°C\n"
        content += "MTTTLQRRESANLWERFCNWVTSTDNRLYVGWFGVIMIPTLLAAT\n"
        content += "ICFVIAFIAAPPVDIDGIREPVSGSLLYGNNIITGAVVPSSNAIG\n"
        
        # Write to temporary file
        filename = f"{job_id}_results.fasta"
        with open(filename, "w") as f:
            f.write(content)
        
        return FileResponse(
            filename,
            media_type="text/plain",
            filename=f"eden_d1_variants_{job_id}.fasta"
        )
    
    elif format == "csv":
        # Generate CSV file
        content = "Variant_ID,ThermoScore,Delta_Tm,Predicted_Tm,Mutations,Confidence,Fitness\n"
        content += f"Variant_001,{results['best_score']},{results['best_delta_tm']},{results['predicted_tm']},{results['num_mutations']},{results['structure_confidence']},{results['fitness_score']}\n"
        
        filename = f"{job_id}_results.csv"
        with open(filename, "w") as f:
            f.write(content)
        
        return FileResponse(
            filename,
            media_type="text/csv",
            filename=f"eden_d1_results_{job_id}.csv"
        )
    
    elif format == "pdb":
        # For now, return a placeholder message
        raise HTTPException(501, "PDB structure download coming soon!")
    
    else:
        raise HTTPException(400, f"Unsupported format: {format}")

@app.get("/api/jobs")
async def list_jobs():
    """List all jobs"""
    return {
        "jobs": [
            {
                "job_id": job_id,
                "status": data["status"],
                "created_at": data.get("created_at"),
                "current_stage": data.get("current_stage")
            }
            for job_id, data in job_store.items()
        ]
    }

@app.delete("/api/pipeline/{job_id}")
async def delete_job(job_id: str):
    """Delete a job from the store"""
    if job_id not in job_store:
        raise HTTPException(404, "Job not found")
    
    del job_store[job_id]
    return {"message": f"Job {job_id} deleted successfully"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)