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
import numpy as np
import os

# Import the real pipeline
from pipeline.master_pipeline_v3 import MasterPipelineV3
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
    n_variants: int = 15  # Default to 15 for initial testing

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

async def run_pipeline_task(job_id: str, sequences: List[str], crop_type: str, n_variants: int):
    """Background task to run the REAL production pipeline"""
    print(f"🔄 Starting PRODUCTION pipeline task for job: {job_id}")
    
    try:
        if job_store[job_id]["status"] == "stopped":
            return
            
        job_store[job_id]["status"] = "running"
        
        # Initialize real pipeline
        config = Config()
        pipeline = MasterPipelineV3(config)
        
        # Progress callback
        async def update_progress(progress: float, message: str):
            if job_store[job_id]["status"] == "stopped":
                raise Exception("Pipeline stopped by user")
            
            job_store[job_id]["progress"] = progress
            job_store[job_id]["current_stage"] = message
            
            # Map progress to stage numbers (1-5)
            if progress < 0.2:
                stage = 1
            elif progress < 0.4:
                stage = 2
            elif progress < 0.6:
                stage = 3
            elif progress < 0.8:
                stage = 4
            else:
                stage = 5
            
            job_store[job_id]["stage_number"] = stage
            print(f"   Stage {stage}/5: {message} ({progress*100:.0f}%)")
        
        # Run REAL pipeline with API calls
        print(f"   Running pipeline with {n_variants} variants targeting {crop_type}")
        
        results = await pipeline.run(
            wt_sequence=sequences[0],
            target_temp=50,  # Target 50°C
            n_variants=n_variants,
            validate_all=False,  # Validate top candidates only for cost
            progress_callback=update_progress
        )
        
        # Process real results
        if results.get('final_variants'):
            best_variant = results['final_variants'][0]
            validation = best_variant.get('validation', {})
            
            # Format mutations for display
            mutations_display = []
            for mut in best_variant.get('mutations', []):
                # Parse mutation string (e.g., "A123V")
                if len(mut) >= 3:
                    from_aa = mut[0]
                    to_aa = mut[-1]
                    pos = mut[1:-1]
                    mutations_display.append({
                        "position": int(pos) if pos.isdigit() else 0,
                        "from": from_aa,
                        "to": to_aa,
                        "effect": "Thermostability enhancement"
                    })
            
            job_store[job_id]["results"] = {
                "variants": len(results['final_variants']),
                "best_score": round(best_variant.get('consensus_score', 0), 2),
                "best_delta_tm": round(validation.get('delta_tm', 0), 1),
                "predicted_tm": round(validation.get('predicted_tm', 0), 1),
                "baseline_tm": round(validation.get('predicted_tm', 0) - validation.get('delta_tm', 0), 1),
                "num_mutations": len(best_variant.get('mutations', [])),
                "structure_confidence": round(validation.get('plddt', 0) / 100.0, 2),
                "fitness_score": round(validation.get('composite_score', 0) / 100.0, 2),
                "mutations": mutations_display,
                "crop_type": crop_type,
                "input_sequence_length": len(sequences[0]),
                "risk_level": best_variant.get('risk_level', 'UNKNOWN'),
                "output_files": results.get('output_files', {}),
                "api_calls": results.get('api_usage', {}).get('total_calls', 0),
                "all_variants": results['final_variants']  # Store all for download
            }
            
            job_store[job_id]["status"] = "complete"
            job_store[job_id]["current_stage"] = "Complete"
            job_store[job_id]["completed_at"] = datetime.now().isoformat()
            
            print(f"✅ Job {job_id} completed successfully")
            print(f"   - Final variants: {len(results['final_variants'])}")
            print(f"   - Best ΔTm: +{validation.get('delta_tm', 0):.1f}°C")
            print(f"   - Risk level: {best_variant.get('risk_level')}")
            print(f"   - API calls made: {results.get('api_usage', {}).get('total_calls', 0)}")
            
        else:
            raise Exception("No variants passed validation criteria")
            
    except Exception as e:
        print(f"❌ Job {job_id} failed with error: {str(e)}")
        import traceback
        traceback.print_exc()
        
        job_store[job_id]["status"] = "failed"
        job_store[job_id]["error"] = str(e)
        job_store[job_id]["current_stage"] = "Failed"
        job_store[job_id]["failed_at"] = datetime.now().isoformat()

@app.post("/api/pipeline/stop/{job_id}")
async def stop_pipeline(job_id: str):
    """Stop a running pipeline job"""
    print(f"🛑 Stop request for job: {job_id}")
    
    if job_id not in job_store:
        raise HTTPException(404, "Job not found")
    
    job_store[job_id]["status"] = "stopped"
    return {"status": "stopped", "job_id": job_id}

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
    
    # Create temp directory if it doesn't exist
    os.makedirs("temp", exist_ok=True)
    
    if format == "fasta":
        # Generate REAL FASTA file with all variants
        filename = f"temp/{job_id}_results.fasta"
        with open(filename, "w") as f:
            for i, variant in enumerate(results.get('all_variants', [])):
                validation = variant.get('validation', {})
                f.write(f">{variant.get('id', f'VAR_{i+1:04d}')} ")
                f.write(f"score={variant.get('consensus_score', 0):.2f} ")
                f.write(f"delta_tm={validation.get('delta_tm', 0):+.1f} ")
                f.write(f"risk={variant.get('risk_level', 'UNKNOWN')}\n")
                f.write(f"{variant['sequence']}\n")
        
        return FileResponse(
            filename,
            media_type="text/plain",
            filename=f"eden_d1_variants_{job_id}.fasta"
        )
    
    elif format == "csv":
        # Generate CSV with detailed results
        filename = f"temp/{job_id}_results.csv"
        with open(filename, "w") as f:
            f.write("Variant_ID,Track,Method,Consensus_Score,Delta_Tm,Predicted_Tm,Risk_Level,Num_Mutations\n")
            
            for i, variant in enumerate(results.get('all_variants', [])):
                validation = variant.get('validation', {})
                f.write(f"{variant.get('id', f'VAR_{i+1:04d')},")
                f.write(f"{variant.get('track', 'unknown')},")
                f.write(f"{variant.get('method', 'unknown')},")
                f.write(f"{variant.get('consensus_score', 0):.3f},")
                f.write(f"{validation.get('delta_tm', 0):.1f},")
                f.write(f"{validation.get('predicted_tm', 0):.1f},")
                f.write(f"{variant.get('risk_level', 'UNKNOWN')},")
                f.write(f"{len(variant.get('mutations', []))}\n")
        
        return FileResponse(
            filename,
            media_type="text/csv",
            filename=f"eden_d1_results_{job_id}.csv"
        )
    
    else:
        raise HTTPException(400, f"Unsupported format: {format}")

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

@app.get("/api/sequences/reference")
async def get_reference_sequences():
    """Get all reference sequences"""
    return {"sequences": REFERENCE_SEQUENCES}

@app.get("/api/jobs")
async def list_jobs():
    """List all jobs"""
    return {
        "jobs": [
            {
                "job_id": job_id,
                "status": data["status"],
                "created_at": data.get("created_at"),
                "current_stage": data.get("current_stage"),
                "progress": data.get("progress", 0),
                "crop_type": data.get("results", {}).get("crop_type") if "results" in data else None
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

    