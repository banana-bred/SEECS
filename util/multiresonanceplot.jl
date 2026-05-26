using Plots
using DelimitedFiles

function collect_files(args)
  isempty(args) && error("Must supply at least one file")
  files = String[]
  for arg in args
    isdir(arg) && error("Must supply files, not directories.")
    isfile(arg) || error("Argument $arg is not a file.")
    push!(files, arg)
  end
  return files
end

function read_energies(path :: AbstractString)
  data = readdlm(path; comments=true, comment_char='#')
  ncols = size(data, 2)
  if ncols == 3
    error("Data at $path is formatted as bound states (ncols = $ncols). No plot to make.")
  elseif ncols != 5
    error("Data at $path is formatted in some unknown way. Was expecting 5 columns, got $ncols")
  end
  reals = Float64.(data[:,2])
  imags = Float64.(data[:,3])
  return reals, imags
end

pathlabel(path) = splitext(basename(path))[1]

function plot_energies(
      files
    ; outfile = "energies.png"
    , kwargs...
  )
  files = collect_files(files)
  isempty(files) && error("No files to plot.")
  default(
    ; xlabel = "Re(E)"
    , ylabel = "Im(E)"
    , title = "SEECS complex energy spectrum"
    , legend = :best
    , aspect_ratio = :equal
    , framestyle = :box
    , markerstrokewidth = 0.3
    , markersize = 4
    , size = (1000, 700)
    , alpha = 0.8
    , dpi = 300
  )
  p=plot()
  for f in files
    r,i = read_energies(f)
    scatter!(p, r, i; label=pathlabel(f), kwargs...)
  end
  hline!(p, [0.0];  label=false, linestyle=:dash,  alpha=0.7, color=:black)
  outfile = "energies.png"
  savefig(p, outfile)
  println("Wrote plot to $outfile")
  return p
end

function main()
  plot_energies(ARGS)
end

if abspath(PROGRAM_FILE) == @__FILE__
  main()
end
