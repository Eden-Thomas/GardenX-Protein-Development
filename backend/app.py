from fastapi import FastAPI, HTTPException, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse  # ADD THIS
from pydantic import BaseModel
from typing import List, Dict, Optional
import json
import asyncio
from datetime import datetime

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

# ADD THIS ROUTE TO SERVE YOUR HTML
@app.get("/")
async def read_root():
    return FileResponse("index.html")

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
    print(f"   - Sequence length: {len(request.sequences[0]) if request.sequences else 0}")
    
    # Validate sequences
    for seq in request.sequences:
        if not all(aa in "ACDEFGHIKLMNPQRSTVWY" for aa in seq.upper()):
            print(f"❌ Invalid amino acid detected in sequence")
            raise HTTPException(400, "Invalid amino acid sequence")
    
    print(f"✅ Sequence validation passed")
    
    # Create job ID
    job_id = f"job_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    print(f"📋 Created job ID: {job_id}")
    
    # Initialize job status
    job_store[job_id] = {
        "status": "queued",
        "progress": 0,
        "created_at": datetime.now().isoformat()
    }
    print(f"💾 Job stored in memory")
    
    # Run pipeline in background
    background_tasks.add_task(
        run_pipeline_task,
        job_id,
        request.sequences,
        request.crop_type,
        request.n_variants
    )
    print(f"🚀 Background task started")
    
    response = PipelineResponse(
        job_id=job_id,
        status="queued",
        message=f"Pipeline started with {len(request.sequences)} sequences"
    )
    print(f"📤 Sending response: {response}")
    
    return response

async def run_pipeline_task(job_id: str, sequences: List[str], crop_type: str, n_variants: int):
    """Background task to run the pipeline"""
    print(f"🔄 Starting pipeline task for job: {job_id}")
    
    try:
        # Check if job was stopped before starting
        if job_store[job_id]["status"] == "stopped":
            print(f"⚠️ Job {job_id} was stopped before starting")
            return
            
        # Update status
        job_store[job_id]["status"] = "running"
        print(f"📊 Job {job_id} status updated to: running")
        
        # Initialize pipeline
        print(f"🔧 Initializing pipeline with Config:")
        print(f"   - Device: {Config.DEVICE if hasattr(Config, 'DEVICE') else 'cpu'}")
        print(f"   - GPU Available: {Config.USE_LOCAL_GPU}")
        print(f"   - API Key Present: {bool(Config.NEUROSNAP_API_KEY)}")
        
        pipeline = D1EngineeringPipeline(Config())
        
        # Run pipeline with progress updates
        print(f"🧬 Running pipeline for {crop_type} with {n_variants} variants")
        results = await pipeline.run_pipeline(
            sequences=sequences,
            crop_type=crop_type,
            n_variants=n_variants,
            progress_callback=lambda p: update_progress(job_id, p)
        )
        
        # Check if job was stopped during execution
        if job_store[job_id]["status"] == "stopped":
            print(f"⚠️ Job {job_id} was stopped during execution")
            return
            
        # Store results
        job_store[job_id]["status"] = "complete"
        job_store[job_id]["results"] = results
        job_store[job_id]["completed_at"] = datetime.now().isoformat()
        print(f"✅ Job {job_id} completed successfully")
        
    except Exception as e:
        print(f"❌ Job {job_id} failed with error: {str(e)}")
        import traceback
        traceback.print_exc()
        job_store[job_id]["status"] = "failed"
        job_store[job_id]["error"] = str(e)

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

@app.get("/api/pipeline/status/{job_id}")
async def get_status(job_id: str):
    """Get the status of a pipeline job"""
    if job_id not in job_store:
        raise HTTPException(404, "Job not found")
    
    return job_store[job_id]