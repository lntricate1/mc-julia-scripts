using NBT, Litematica, MinecraftDataStructures

"""
    setup(max_power::Float32)

Precalculates the relationships between rays and blocks for explosions of at most max_power.

# Returns
- `points_per_block::Vector{Int}`: The number of rays affecting each block.
- `ray_indices::Vector{Tuple{Int, Float32, Float32}}`: For each block, `n` entries (given by `points_per_block`) of `(ray_index, dampening_air, dampening_block)`.
- `block_counts::Vector{Int}`: The amount of each block after removing symmetry.
- `dict::Dict{Int, Tuple{Int, Int, Int}}`: Maps block index (as in `points_per_block` and `block_counts`) to block position.
"""
function setup(max_power::Float32)
  rays = Vector{Tuple{Float64, Float64, Float64}}(undef, 36)
  i = 0
  for y in 8:15, x in 8:y
    X, Y, Z = (x, y, 15) ./ 15 .* 2 .- 1
    mag = sqrt(X^2 + Y^2 + Z^2)
    rays[i+=1] = (X, Y, Z) ./ mag .* 0.3f0
  end

  max_strength = 1.3f0max_power
  ray_indices = Vector{Pair{Tuple{Int, Int, Int}, Tuple{Int, Float32, Float32}}}[]
  points_per_block = Dict{Tuple{Int, Int, Int}, Int}()
  for (ind, ray) in enumerate(rays)
    pos = (1.5, 1.5, 1.5)
    strength = max_strength
    dampening_air = 0f0
    dampening_block = 0f0
    depth = 1
    prev_depth = 1
    prevx = 1
    prevy = 1
    prevz = 1
    while strength > 0.09f0
      x = trunc(Int, pos[1])
      y = trunc(Int, pos[2])
      z = trunc(Int, pos[3])
      if prevx != x || prevy != y || prevz != z
        while length(ray_indices) < depth
          push!(ray_indices, Pair{Tuple{Int, Int, Int}, Tuple{Int, Float32, Float32}}[])
        end
        push!(ray_indices[prev_depth], (prevx, prevy, prevz) => (ind, dampening_air, dampening_block))
        if haskey(points_per_block, (prevx,prevy,prevz))
          points_per_block[(prevx,prevy,prevz)] += 1
        else
          points_per_block[(prevx,prevy,prevz)] = 1
        end
        dampening_air = 0f0
        dampening_block = 0f0
        prev_depth = depth
        prevx = x
        prevy = y
        prevz = z
      end
      pos = pos .+ ray
      strength -= 0.22500001f0
      dampening_block += 0.09f0
      dampening_air += 0.22500001f0
      dampening_block += 0.22500001f0
      depth += 1
    end
  end

  out_ray_indices = Vector{Tuple{Int, Float32, Float32}}(undef, sum(values(points_per_block)))
  out_points_per_block = Vector{Int}(undef, length(points_per_block))
  out_block_counts = Vector{Int}(undef, length(points_per_block))
  out_dict = Dict(i => k for (i, k) in enumerate(keys(points_per_block)))
  block_indices = Dict{Tuple{Int, Int, Int}, Int}()
  next_index = 1
  j = 0
  for datas in ray_indices
    for pair in datas
      index = pair.first
      data = pair.second
      if haskey(block_indices, index)
        i = block_indices[index] += 1
      else
        n = points_per_block[index]
        i = next_index
        block_indices[index] = next_index
        next_index += n
        out_points_per_block[j+=1] = n
        x, y, z = index
        out_block_counts[j] =
          x == y == z == 1 ? 1 :
          x == y == 1 ? 6 :
          x == y == z ? 8 :
          y == z && x == 1 ? 12 :
          x <= y <= z ? x == y || y == z || x == 1 ? 24 : 48 : 0
        out_dict[j] = index
      end
      out_ray_indices[i] = data
    end
  end
  return out_points_per_block, out_ray_indices, out_block_counts, out_dict
end

"""
    gen_points(blocks::UInt64, ray_indices::Vector{Tuple{Int, Float32, Float32}}, points_per_block::Vector{Int})

Generates ray strengths at key points for a given block layout `blocks`.
"""
function gen_points(blocks::UInt64, ray_indices::Vector{Tuple{Int, Float32, Float32}}, points_per_block::Vector{Int})
  rays = zeros(Float32, 36)
  points = Vector{Float32}(undef, length(ray_indices))
  b = 0x0000000000000001
  i = 0
  for p in points_per_block
    if blocks & b != 0x0000000000000000
      for _ in 1:p
        rayind, damp = ray_indices[i+=1]
        points[i] = -1e10
        rays[rayind] -= damp
      end
    else
      for _ in 1:p
        rayind, _, damp = ray_indices[i+=1]
        points[i] = rays[rayind] - 0.09f0
        rays[rayind] -= damp
      end
    end
    b *= 0x0000000000000002
  end
  return points
end

cum(x::Float32) = 2f0/9f0 * (
  x <= 4f0*0.7f0 ? 0f0 :
  x <= 4f0*1.3f0 ? x*(log(x) - log(4f0*0.7f0) - 1) + 4f0*0.7f0 :
  x <= 11.5f0*0.7f0 ? x*log(1.3f0/0.7f0) + 4f0(0.7f0-1.3f0) :
  x <= 11.5f0*1.3f0 ? x*(log(11.5f0*1.3f0) - log(x) + 1) + 4f0(0.7f0-1.3f0) - 11.5f0*0.7f0 :
  9f0/2f0
)

"""
    X_cart(blocks::UInt64, counts::Vector{Float32}, strength_mags_inv::Vector{Float32}, rays::Vector{Float32}, ray_indices::Vector{Tuple{Int, Float32, Float32}}, points_per_block::Vector{Int}, block_counts::Vector{Int})

Calculates average blocks broken by a cart for the block layout `blocks`. See [`bruteforce_chamber_cart`](@ref) for an example of how to use this.
"""
function X_cart(blocks::UInt64, counts::Vector{Float32}, strength_mags_inv::Vector{Float32}, rays::Vector{Float32}, ray_indices::Vector{Tuple{Int, Float32, Float32}}, points_per_block::Vector{Int}, block_counts::Vector{Int})
  rays .= 0f0
  b = 0x0000000000000001
  i = 0
  K = 1f0 - 1.3f0/0.6f0
  total = 0f0
  total1 = 0f0
  for (j, p) in enumerate(points_per_block)
    if blocks & b != 0x0000000000000000
      for _ in 1:p
        rayind, damp = ray_indices[i+=1]
        rays[rayind] -= damp
      end
    else
      counts .= 1f0
      counts1 = 1f0
      for _ in 1:p
        rayind, _, damp = ray_indices[i+=1]
        s = rays[rayind] - 0.09f0
        # println(clamp.(K .- s.*strength_mags_inv, 0f0, 1f0), ", ", s, ", ", cum(11.5f0*1.7f0 + s))
        # counts1 *= cum(11.5f0*1.7f0 + s)
        counts1 *= cum(-s)
        counts .*= clamp.(K .- s.*strength_mags_inv, 0f0, 1f0)
        rays[rayind] -= damp
      end
      counts .-= 1f0
      total += sum(counts)*block_counts[j]
      total1 += (1-counts1)*block_counts[j]
      a = -sum(counts)/length(counts)
      c = 1-counts1
      println(a-c, ", ", a, ", ", c)
    end
    b *= 0x0000000000000002
  end
  println((-total/length(counts), total1))
  return -total/length(counts)
end

"""
    X(blocks::UInt64, counts::Vector{Float32}, strength_mags_inv::Vector{Float32}, rays::Vector{Float32}, ray_indices::Vector{Tuple{Int, Float32, Float32}}, points_per_block::Vector{Int}, block_counts::Vector{Int})

Calculates average blocks broken for the block layout `blocks`. See [`bruteforce_chamber`](@ref) for an example of how to use this.
"""
function X(blocks::UInt64, strength_mag_inv::Float32, rays::Vector{Float32}, ray_indices::Vector{Tuple{Int, Float32, Float32}}, points_per_block::Vector{Int}, block_counts::Vector{Int})
  rays .= 0f0
  b = 0x0000000000000001
  i = 0
  K = 1f0 - 1.3f0/0.6f0
  total = 0f0
  for (j, p) in enumerate(points_per_block)
    if blocks & b != 0x0000000000000000
      for _ in 1:p
        rayind, damp = ray_indices[i+=1]
        rays[rayind] -= damp
      end
    else
      count = 1f0
      for _ in 1:p
        rayind, _, damp = ray_indices[i+=1]
        s = rays[rayind] - 0.09f0
        count *= clamp(K - s*strength_mag_inv, 0f0, 1f0)
        rays[rayind] -= damp
      end
      total += (1f0-count)*block_counts[j]
    end
    b *= 0x0000000000000002
  end
  return total
end

"""
    bruteforce_chamber_cart(iterations::Int=200, iterator::StepRange{UInt, UInt}=UInt(1):UInt(2):UInt(2^17), mask::UInt=UInt(1))

Finds the optimal minecart blast chamber layout.

# Arguments
- `mask::UInt`: the blocks that are forced to be air.
- `iterator::StepRange{UInt, UInt}`: specifies how to iterate over block layouts, as `UInt`s where 1 bits are air.
- `iterations::Int`: how many explosion powers to use for averaging.

# Examples
```julia-repl
julia> @time bruteforce_chamber_cart()
  7.995023 seconds (248 allocations: 269.500 KiB)
(3209.7358f0, 0x000000000000fdff)

julia> @time bruteforce_chamber_cart(200, UInt(0):UInt(2^16):UInt(2^32), 0x000000000000fdff) # After a basic test you can mask out unnecessary bits to look further
  7.708847 seconds (252 allocations: 272.641 KiB)
(3225.819f0, 0x000000000458fdff)
```
"""
function bruteforce_chamber_cart(iterations::Int=200, iterator::StepRange{UInt, UInt}=UInt(1):UInt(2):UInt(2^17), mask::UInt=UInt(1))
  points_per_block, ray_indices, block_counts, dict = setup(11.5f0)
  strength_mags_inv = inv.(0.6f0range(4f0, 11.5f0, iterations))
  # Dummy stuff
  counts = zeros(Float32, iterations)
  rays = zeros(Float32, 36)

  best = 0f0
  bestblocks = 0x0000000000000000
  for blocks in iterator
    count = X_cart(blocks | mask, counts, strength_mags_inv, rays, ray_indices, points_per_block, block_counts)
    if count > best
      best = count
      bestblocks = blocks | mask
    end
  end
  return best, bestblocks
end

"""
    bruteforce_chamber(;iterations::Int=200, iterator::StepRange{UInt, UInt}=UInt(1):UInt(2):UInt(2^17), mask::UInt=UInt(1))

Finds the optimal minecart blast chamber layout.

# Arguments
- `mask::UInt`: the blocks that are forced to be air.
- `iterator::StepRange{UInt, UInt}`: specifies how to iterate over block layouts, as `UInt`s where 1 bits are air.
- `iterations::Int`: how many explosion powers to use for averaging.

# Examples
```julia-repl
julia> @time bruteforce_chamber()
  0.257991 seconds (101 allocations: 66.500 KiB)
(604.5703f0, 0x00000000000000ff)

julia> @time bruteforce_chamber(11.5f0, UInt(1):UInt(2):UInt(2^19))
  1.369937 seconds (244 allocations: 267.688 KiB)
(7076.188f0, 0x000000000007fdff)
```
"""
function bruteforce_chamber(power::Float32=4f0, iterator::StepRange{UInt, UInt}=UInt(1):UInt(2):UInt(2^19), mask::UInt=UInt(1))
  points_per_block, ray_indices, block_counts, dict = setup(power)
  strength_mag_inv = inv(0.6f0power)
  # Dummy stuff
  rays = zeros(Float32, 36)

  best = 0f0
  bestblocks = 0x0000000000000000
  for blocks in iterator
    count = X(blocks | mask, strength_mag_inv, rays, ray_indices, points_per_block, block_counts)
    if count > best
      best = count
      bestblocks = blocks | mask
    end
  end
  return best, bestblocks
end

"""
    show_solution(blocks::UInt, power::Float32=4f0)

Displays a 1/48th symmetrical slice of the blast chamber efficiency using a Makie `voxels` plot.
"""
function show_solution(blocks::UInt, power::Float32=4f0, colormap=:inferno)
  points_per_block, ray_indices, block_counts, dict = setup(power)
  points = gen_points(blocks, ray_indices, points_per_block)

  maxpos = (1,1,1)
  for b in values(dict)
    maxpos = max.(maxpos, b)
  end
  strength_mag = 0.6f0power
  K = 1f0 - 1.3f0/0.6f0
  count = 1f0
  i = 1
  j = 1
  total = 0f0
  data = zeros(Float32, maxpos)
  for p in points
    count *= clamp(K - p/strength_mag, 0f0, 1f0)
    if (j += 1) > points_per_block[i]
      total += block_counts[i] * (1f0 - count)
      data[dict[i]...] = 1f0 - count
      count = 1f0
      i += 1
      j = 1
    end
  end
  println(total)
  voxels(data, is_air=(==(0f0)), colormap=colormap)
end

function _makethedata!(power::Float32, blocks::BitArray{3}, block_counts::Array{Int, 3}, rays::Vector{Tuple{Float64, Float64, Float64}}, pos::Tuple{Float64, Float64, Float64}, data::Array{Float32, 3}, raycounts::Vector{Int})
  strength_mag_inv = inv(0.6f0power)
  data .= 1f0
  K = 1f0 - 1.3f0/0.6f0
  for (i, r) in enumerate(rays)
    strength = 1.3f0power
    damp = 0f0
    pos_ = pos
    px = -1; py = -1; pz = -1
    raycount = 0
    while strength > 0.09f0
      x, y, z = trunc.(Int, pos_)
      raycount += 1
      if blocks[x,y,z]
        strength -= 0.09f0
        damp += 0.09f0
        if x != px || y != py || z != pz
          # print(K + strength*strength_mag_inv, ", ")
          # print(1.3f0/0.6f0 - strength*strength_mag_inv, ", ")
          # print(K + damp*strength_mag_inv, ", ")
          # data.val[x,y,z] *= clamp(1.3f0/0.6f0 - strength*strength_mag_inv, 0f0, 1f0)
          data[x,y,z] *= clamp(K + damp*strength_mag_inv, 0f0, 1f0)
          # data.val[x,z,y] *= clamp(K + damp*strength_mag_inv, 0f0, 1f0)
          # data.val[y,x,z] *= clamp(K + damp*strength_mag_inv, 0f0, 1f0)
          # data.val[y,z,x] *= clamp(K + damp*strength_mag_inv, 0f0, 1f0)
          # data.val[z,x,y] *= clamp(K + damp*strength_mag_inv, 0f0, 1f0)
          # data.val[z,y,x] *= clamp(K + damp*strength_mag_inv, 0f0, 1f0)
        end
      end
      strength -= 0.22500001f0
      damp += 0.22500001f0
      px = x; py = y; pz = z;
      pos_ = pos_ .+ r
    end
    raycounts[i] = raycount
  end
  # println(data.val)
  data .= 1f0 .- data
  return dot(block_counts, data)
end

"""
    draw_wiki_thing(blocks::UInt, power::Float32=4f0)

Draws a diagram for the Minecraft wiki explosion page.
"""
function show_solution2(blocks::BitArray{3}, colormap::Symbol=:inferno, power::Float32=4f0)
  M = size(blocks)[1]
  data = zeros(Float32, M, M, M)
  raycounts = zeros(Int, 169)
  raydirs = Vector{Tuple{Float64, Float64, Float64}}(undef, 169)
  raydirsnorm = Vector{Tuple{Float64, Float64, Float64}}(undef, 169)
  i = 0
  for x in 8:15, y in 8:15, z in 8:15 if x==15 || y==15 || z==15
    X, Y, Z = (x, y, z) ./ 15 .* 2 .- 1
    mag = sqrt(X^2 + Y^2 + Z^2)
    raydirs[i+=1] = (X, Y, Z).*M
    raydirsnorm[i] = (X, Y, Z) ./ mag .* 0.3f0
  end end
  fig = Figure()
  ax = LScene(fig[1,1], show_axis=false)
  # sg = SliderGrid(
  #   fig[0, 1],
  #   (label = "X", range = -0.5:0.025:0.5, startvalue = 0.0),
  #   (label = "Y", range = -0.5:0.025:0.5, startvalue = 0.0),
  #   (label = "Z", range = -0.5:0.025:0.5, startvalue = 0.0),
  # )
  # Pos = lift((x,y,z) -> (x,y,z), sg.sliders[1].value, sg.sliders[2].value, sg.sliders[3].value)
  Pos = Observable((0.0, 0.0, 0.0))
  Raydirs = lift(p -> [r .+ p for r in raydirs], Pos)
  block_counts = [
    x == y == z == 1 ? 1 :
    x == y == 1 || y == z == 1 || x == z == 1 ? 2 :
    x == 1 || y == 1 || z == 1 ? 4 :
    8 for x in 1:M, y in 1:M, z in 1:M]

  Data = lift(Pos) do pos
    _makethedata!(power, blocks, block_counts, raydirsnorm, pos .+ 1.5, data, raycounts)
    return data
  end

  _makethedata!(power, blocks, block_counts, raydirsnorm, (1.5,1.5,1.5), data, raycounts)
  plt = voxels!(ax, -0.5 .. M-0.5, -0.5 .. M-0.5, -0.5 .. M-0.5, Data, colorrange=(0.0,1.0), colormap=colormap, is_air=(==(0f0)))
  Colorbar(fig[2,1], plt, nsteps=256, ticks=0:0.1:1, vertical=false, halign=:center, label="Block break probability")
  Sc1 = lift(pos -> [pos .+ p.*i for (j, p) in enumerate(raydirsnorm) for i in 1:raycounts[j]], Pos)
  sc1 = scatter!(Sc1, overdraw=true, color=:red, markersize=2.0, label="Ray points")
  sc2 = scatter!(Raydirs, overdraw=true, color=:black, markersize=2.0, label="Ray directions")
  for i in eachindex(raydirs)
    lines!(lift((p, r) -> [p.+r[i]./2, r[i]], Pos, Raydirs), color=:black, linewidth=0.2)
  end
  sc3 = scatter!(lift((p) -> [p], Pos), strokecolor=:black, strokewidth=2, color=:magenta, overdraw=true)
  axislegend(ax, [sc1 => (;markersize=20), sc2 =>(;markersize=20), sc3 =>(;markersize=20)], ["Ray points", "Ray directions", "Explosion"])

  cam = cam3d!(ax)
  rotate_cam!(ax.scene, 3pi/16, 0, 0)
  rotate_cam!(ax.scene, 0, 10pi/8, 0)
  # zoom!(ax.scene, cam, 0.3)
  # zoom!(ax.scene, cam, 0.3)
  iterations = 200
  record(fig, "rotating thingy.mp4", range(0, 2pi, iterations)[1:end-1], framerate=10) do rot
    rotate_cam!(ax.scene, 0.0, 2pi/iterations, 0.0)
  end
  # iterations = 200
  # record(fig, "balls.gif", range(0,4pi,iterations), framerate=10) do t
  #   rotate_cam!(ax.scene, 0.0, 2pi/iterations, 0.0)
  #   r = 0.4
  #   pos = (r*cos(t), r*sin(t), r*sin(t))
  #   x, y, z = pos
  #   _makethedata!(power, blocks, block_counts, raydirsnorm, (1.5 + x, 1.5 + y, 1.5 + z), data, raycounts)
  #   Makie.local_update(plt, :, :, :)
  #   Pos[] = pos
  # end
  return fig, ax, cam
end

"""
    show_solution3(blocks::UInt, power::Float32=4f0)

Displays the blast chamber efficiency for all blocks using a Makie `voxels` plot.
"""
function show_solution3(blocks::UInt, power::Float32=4f0, colormap=:inferno)
  points_per_block, ray_indices, block_counts, dict = setup(power)
  points = gen_points(blocks, ray_indices, points_per_block)

  maxpos = (1,1,1)
  for b in values(dict)
    maxpos = max.(maxpos, b)
  end
  strength_mag = 0.6f0power
  K = 1f0 - 1.3f0/0.6f0
  count = 1f0
  i = 1
  j = 1
  total = 0f0
  M = maximum(maxpos)
  data = zeros(Float32, 2M-1, 2M-1, 2M-1)
  for p in points
    count *= clamp(K - p/strength_mag, 0f0, 1f0)
    if (j += 1) > points_per_block[i]
      total += block_counts[i] * (1f0 - count)
      X, Y, Z = dict[i]
      for x in (M-X+1, M+X-1), y in (M-Y+1, M+Y-1), z in (M-Z+1, M+Z-1)
        data[x,y,z] = 1f0 - count
        data[x,z,y] = 1f0 - count
        data[y,x,z] = 1f0 - count
        data[y,z,x] = 1f0 - count
        data[z,x,y] = 1f0 - count
        data[z,y,x] = 1f0 - count
      end
      count = 1f0
      i += 1
      j = 1
    end
  end
  println(total)
  voxels(data, is_air=(==(0f0)), colormap=colormap)
end

"""
    make_schem(blocks::UInt[, metadata::TagCompound])

Makes a `Litematic` of the blast chamber layout.
"""
function make_schem(blocks::UInt, metadata::TagCompound=TagCompound([
  "TimeCreated" => round(Int, 1000time()),
  "TimeModified" => round(Int, 1000time()),
  "RegionCount" => 1,
  "Author" => "funny bruteforcer",
  "Name" => "good(?) blast chamber layout"
  ]))
  dict = setup(11.5f0)[4]
  maxpos = (1,1,1)
  for b in values(dict)
    maxpos = max.(maxpos, b)
  end
  M = maximum(maxpos)
  data = ones(Block, 2M-1, 2M-1, 2M-1)
  b = 0x0000000000000001
  for i in 1:64
    if b & blocks != 0x0000000000000000
      X, Y, Z = dict[i]
      for x in (M-X+1, M+X-1), y in (M-Y+1, M+Y-1), z in (M-Z+1, M+Z-1)
        data[x,y,z] = zero(Block)
        data[x,z,y] = zero(Block)
        data[y,x,z] = zero(Block)
        data[y,z,x] = zero(Block)
        data[z,x,y] = zero(Block)
        data[z,y,x] = zero(Block)
      end
    end
    b *= 0x0000000000000002
  end
  push!(metadata.data, "TotalBlocks" => sum(data .!= zeros(Block, size(data))))
  push!(metadata.data, "EnclosingSize" => TagCompound(["x" => Int32(2M-1), "y" => Int32(2M-1), "z" => Int32(2M-1)]))
  push!(metadata.data, "TotalVolume" => (2M-1)^3)
  return Litematic(2586, metadata, [Region(data)])
end
