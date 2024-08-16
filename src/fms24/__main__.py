import argparse
from .bomd import BOMD
from .bundle import Bundle

def main():
    parser = argparse.ArgumentParser(description="Run BOMD simulation")
    parser.add_argument('--config', type=str, default='config.yaml', help='Path to configuration file')
    args = parser.parse_args()

    # Initialize BOMD
    bomd = BOMD(config=args.config)
    bomd.initialize()

    # Initialize Bundle
    bundle = Bundle(config=args.config)

    # Run simulation
    # This is a placeholder for your actual simulation loop
    for step in range(100):  # Replace with actual number of steps
        bomd.prop_velocity_verlet(dt=0.001)  # Replace with actual time step
        bundle.propagate_tbfs(dt=0.001)
        
        # Add logic for dynamically adding/removing TBFs
        
        if step % 10 == 0:  # Write trajectory every 10 steps
            bomd.write_traj_dump()

    print("Simulation completed")

if __name__ == "__main__":
    main()