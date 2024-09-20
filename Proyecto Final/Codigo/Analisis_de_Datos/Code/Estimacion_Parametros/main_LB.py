import sys
import os
from LB import LatticeBoltzmann, LoadData, LoadDataFires, PrintData, PrintFrame

def main():
    # Parse command-line arguments
    if len(sys.argv) != 2:
        print("Usage: python script.py <diffusion_coefficient>")
        sys.exit(1)

    D = float(sys.argv[1])  # Diffusion coefficient
    iter_per_hour = 5
    t_hour = 240

    # Simulation parameters
    tframe = iter_per_hour
    tmax = iter_per_hour * t_hour
    delta_t = 1

    # Initial conditions
    rho0 = 0.00001
    Ux0 = 0.0
    Uy0 = 0.0

    # Create LatticeBoltzmann instance
    air = LatticeBoltzmann()

    # Load data
    air.LoadData("velocity.txt")
    air.LoadDataFires("coordenadasfuentesgrilla1410.txt")

    # Calculate relaxation time tau
    tau = (D / delta_t * air.Cs2) + 0.5
    Utau = 1.0 / tau
    UmUtau = 1 - Utau

    # Start simulation
    air.Start(rho0, Ux0, Uy0)

    # Main simulation loop
    for t in range(tmax + 1):
        air.Collision()
        air.ImposeFields(t)
        air.Advection()

        # Print results every tframe steps
        if t % tframe == 0:
            # Create data directory if it doesn't exist
            data_dir = "data"
            if not os.path.exists(data_dir):
                os.makedirs(data_dir)

            # Generate filename for density data
            filename = f"data/density_{t:03d}.dat"

            # Save density data
            air.PrintData(filename, t)

            # Generate frame using Gnuplot
            air.PrintFrame(t)

            # Print percentage of completion
            print(f"Porcentaje de avance: {(t * 100) / tmax}%")

if __name__ == "__main__":
    main()
