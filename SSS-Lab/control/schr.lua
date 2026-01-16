-- Define the grid and operator
local domain = Grid.new(256, 256) -- 2D Domain
local H = Schrodinger.new(domain)

-- Configure Potential: Harmonic Oscillator with perturbation
H:set_potential(function(x, y)
    return (x^2 + y^2) + 0.1 * math.random() 
end)

-- Solve for the first 10 eigenvalues
local spectrum = H:compute_spectrum(10)

-- Compute Doubling Index at origin
local d_idx = H:analyze_doubling_index(0, 0, 0.5) -- radius 0.5
print("First eigenvalue:", spectrum[0])
print("Doubling index at origin:", d_idx)