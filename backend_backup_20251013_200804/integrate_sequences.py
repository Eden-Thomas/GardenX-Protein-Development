"""
integrate_sequences.py
Integrate downloaded sequences into analysis_engine.py
"""

from pathlib import Path
import sys

def integrate_new_sequences():
    """Add downloaded sequences to REFERENCE_SEQUENCES"""
    
    print("="*70)
    print("INTEGRATING NEW SEQUENCES INTO PLATFORM")
    print("="*70)
    
    # Check if downloaded sequences exist
    seq_file = Path("data/external_databases/downloaded_sequences.py")
    if not seq_file.exists():
        print("\n✗ Error: downloaded_sequences.py not found")
        print("Please run: python backend/download_datasets.py first")
        return False
    
    # Import the downloaded sequences
    sys.path.insert(0, str(Path("data/external_databases")))
    from downloaded_sequences import ADDITIONAL_SEQUENCES
    
    print(f"\n✓ Found {len(ADDITIONAL_SEQUENCES)} new sequences")
    
    # Read current analysis_engine.py
    engine_file = Path("backend/analysis_engine.py")
    if not engine_file.exists():
        print("✗ Error: backend/analysis_engine.py not found")
        return False
    
    with open(engine_file, 'r') as f:
        content = f.read()
    
    # Check if already integrated
    if 'ADDITIONAL_SEQUENCES' in content:
        print("\n✓ Sequences already integrated!")
        return True
    
    # Find REFERENCE_SEQUENCES dictionary
    if 'REFERENCE_SEQUENCES = {' not in content:
        print("✗ Error: Could not find REFERENCE_SEQUENCES in analysis_engine.py")
        return False
    
    # Add import at top
    import_line = "\nfrom pathlib import Path\nimport sys\n# Import additional sequences\ntry:\n    sys.path.insert(0, str(Path('data/external_databases')))\n    from downloaded_sequences import ADDITIONAL_SEQUENCES\nexcept ImportError:\n    ADDITIONAL_SEQUENCES = {}\n"
    
    # Find imports section and add
    if "import numpy as np" in content:
        content = content.replace(
            "import numpy as np",
            f"import numpy as np{import_line}"
        )
    
    # Add to REFERENCE_SEQUENCES
    # Find the closing brace of REFERENCE_SEQUENCES
    ref_seq_end = content.find("}", content.find("REFERENCE_SEQUENCES"))
    
    if ref_seq_end != -1:
        # Add merge line after REFERENCE_SEQUENCES definition
        merge_code = "\n\n# Merge with additional downloaded sequences\nREFERENCE_SEQUENCES.update(ADDITIONAL_SEQUENCES)\n"
        content = content[:ref_seq_end + 1] + merge_code + content[ref_seq_end + 1:]
    
    # Write back
    with open(engine_file, 'w') as f:
        f.write(content)
    
    print(f"\n✓ Integration complete!")
    print(f"  - Added import for ADDITIONAL_SEQUENCES")
    print(f"  - Merged {len(ADDITIONAL_SEQUENCES)} new sequences")
    print(f"  - Total sequences now available: {5 + len(ADDITIONAL_SEQUENCES)} (5 original + {len(ADDITIONAL_SEQUENCES)} new)")
    
    return True


def verify_integration():
    """Verify integration worked"""
    
    print("\n" + "="*70)
    print("VERIFYING INTEGRATION")
    print("="*70)
    
    try:
        # Import and check
        from analysis_engine import REFERENCE_SEQUENCES
        
        print(f"\n✓ Successfully imported REFERENCE_SEQUENCES")
        print(f"✓ Total sequences available: {len(REFERENCE_SEQUENCES)}")
        print(f"\nAvailable organisms:")
        for org in sorted(REFERENCE_SEQUENCES.keys()):
            seq_len = len(REFERENCE_SEQUENCES[org])
            print(f"  - {org}: {seq_len} amino acids")
        
        return True
        
    except Exception as e:
        print(f"\n✗ Verification failed: {e}")
        return False


def main():
    """Run integration"""
    
    success = integrate_new_sequences()
    
    if success:
        verify_integration()
        
        print("\n" + "="*70)
        print("READY TO RUN ANALYSIS!")
        print("="*70)
        print("\nNext step: python backend/master_pipeline_v1.py")
        print("="*70)
    else:
        print("\n✗ Integration failed. Please check errors above.")


if __name__ == "__main__":
    main()