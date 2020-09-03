# Convergence-divergence game.
# 
# Henri Kauhanen 2020


# Dependencies
using DelimitedFiles
using DifferentialEquations
using LinearAlgebra
using Random


# Time derivative; u is x in the paper, p holds the parameters, and
# t is time. Modifies du in place.
function verger!(du, u, p, t)
  # this makes for easier reading
  S = p.S
  tS = p.tS
  m = p.m
  tm = p.tm
  x = u
  tx = 1 .- u

  # social alignment matrices with zero diagonals
  S0 = S
  S0[diagind(S0)] .= 0
  tS0 = tS
  tS0[diagind(tS0)] .= 0

  # fitnesses
  f = x .* (diag(S) .+ S0*tx)
  tf = tx .* (diag(tS) .+ tS0*x)

  # time derivative
  du .= (tx .- tm).*x.*f .- (x .- m).*tx.*tf
end


# Parameter sweep for asymmetric three-population case
function popthreesweep(outfile)
  # open outfile and write CSV header
  open(outfile, "w") do io
    write(io, "sigma,upsilon,k,x1,x2,x3\n")

    # loop through grid
    for sigma = 0.2:0.2:5
      for upsilon in [0.2, 1, 5]
        for k = 0.2:0.2:5
          S = [1 sigma (1/k)*sigma;
               sigma 1 k*sigma;
               upsilon upsilon 1]
          for iter = 1:100
            result = popthree(S)
            result = result(20.0)
            output = [sigma, upsilon, k, result[1], result[2], result[3]]'
            writedlm(io, output, ",")
          end
        end
      end
    end
  end
end


# Solutions for one parameter combination for asymm. 3-pop. case.
function popthreesweepone(outfile)
  # open outfile and write CSV header
  open(outfile, "w") do io
    write(io, "sigma,upsilon,k,id,t,x1,x2,x3\n")

    sigma = 4
    upsilon = 1
    k = 4
    S = [1 sigma (1/k)*sigma;
         sigma 1 k*sigma;
         upsilon upsilon 1]

    for id = 1:100
      result = popthree(S)
      #for t in 0.0:0.1:20.0
      for t in (10^y for y in range(log10(0.1), log10(20.0), length=100))
        output = [sigma, upsilon, k, id, t, result(t)[1], result(t)[2], result(t)[3]]'
        writedlm(io, output, ",")
      end
    end
  end
end


# Numerical solution of one particular instance of the 3-population game
# (jocks, burnouts and adults), returning the solution
function popthree(S)
  # parameters
  m = 0.13
  tm = 0.23

  # parameters collected in a tuple
  p = (S=S, tS=S, m=m, tm=tm)

  starter = rand(3)
  starter[3] = (1/4)*rand()

  maxtime = 20.0
  prob = ODEProblem(verger!, starter, (0.0, maxtime), p)
  solve(prob, DP5())
end

