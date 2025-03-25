abstract type Pearl end
# abstract type Arrow end
abstract type PearlOld end

# function tick_y(::Type{Arrow}, p::Float64, v::Float64)
#   p += v
#   v *= 0.99f0
#   v -= 0.05f0
#   return p, v
# end

function tick_y(::Type{Pearl}, p::Float64, v::Float64)
  v -= 0.03
  v *= 0.99f0
  p += v
  return p, v
end

function tick_y(::Type{PearlOld}, p::Float64, v::Float64)
  p += v
  v *= 0.99f0
  v -= 0.03f0
  return p, v
end

@inline function _nearestblockpos(pos::Float64)
  block = round(16pos)/16
  return block, pos - block
end

@inline function _nearestblockpos(::Type{Pearl}, pos::Float64, explosionheight::Float64)::Tuple{Float64, Float64}
  return _nearestblockpos(pos + 0.2125f0 - explosionheight)
end

@inline function _nearestblockpos(::Type{PearlOld}, pos::Float64, explosionheight::Float64)::Tuple{Float64, Float64}
  return _nearestblockpos(pos + 0.2125f0 - explosionheight)
end

# @inline function _nearestblockpos(::Type{Arrow}, pos::Float64, explosionheight::Float64)::Tuple{Float64, Float64}
#   return _nearestblockpos(pos + 0.13f0 - explosionheight)
# end

struct BounceIndices{N}
  indices::NTuple{N, Int}
  pos::Float64
  vel::Float64
  entity_type::Type
end

struct BounceIndex{N}
  I::NTuple{N, Int}
  P::NTuple{N, Float64}
  V::NTuple{N, Float64}
end

@inline function __inc(state::NTuple{N, Int}, pos::NTuple{N, Float64}, vel::NTuple{N, Float64}, indices::NTuple{N, Int}, entity_type::Type) where N
  I = first(state)
  ts, tp, tv = Base.tail(state), Base.tail(pos), Base.tail(vel)
  if I < first(indices)
    NP, NV = tick_y(entity_type, first(pos), first(vel))
    return true, (I + 1, ts...), (NP, tp...), (NV, tv...)
  end
  first_zero, I, P, V = __inc(ts, tp, tv, Base.tail(indices), entity_type)
  FP = first(P)
  NP = first_zero ? FP + _slime_bounce_pos(entity_type, FP) : FP
  return false, (1, I...), (NP, P...), (1., V...)
end

@inline function Base.iterate(iter::BounceIndices{N}) where N
  BI = BounceIndex(ntuple(i -> 1, Val(N)), ntuple(i -> iter.pos, Val(N)), ntuple(i -> iter.vel, Val(N)))
  return BI, BI
end

@inline function Base.iterate(iter::BounceIndices{N}, state::BounceIndex{N}) where N
  state.I == iter.indices && return nothing
  _, I, P, V = __inc(state.I, state.P, state.V, iter.indices, iter.entity_type)
  next = BounceIndex(I, P, V)
  return next, next
end

@inline Base.eltype(::Type{BounceIndices{N}}) where N = BounceIndex{N}
@inline Base.length(iter::BounceIndices{N}) where N = prod(iter.indices)

"""
    get_slime_bounces(pos, vel, ticks, threshold; entity_type, explosion_height, min_pos, max_pos, max_ticks)

Simulates all possible slime bounces within the limits specified, and returns only those which are closer to the perfect alignment than `threshold`.

Returns (pos::Vector{Float64}, vel::Vector{Float64}, ticks::Vector{Tuple}, delta::Vector{Float64}, block::Vector{Float64}).

See also [`sim_bounces`](@ref).

# Arguments
- `pos::Float64`: The starting Y position.
- `vel::Float64`: The starting Y velocity.
- `ticks::Tuple`: A list of maximum bounce lengths, in ticks, written (bounce 1 length, bounce 2 length, ..., bounce N length).
- `threshold::Float64`: The biggest acceptable distance from the perfect alignment.
- `entity_type::Type`: The entity type to bounce. Supported types: `Pearl`, `PearlOld`.
- `explosion_height::Float64=0.061250001192092896`: The explosion height of the explosive.
- `min_pos::Float64=pos`: The minimum projectile arrival position.
- `max_pos::Float64=256.0`: The minimum projectile arrival position.
- `max_ticks::Int=sum(ticks)`: The maximum total ticks.

# Examples
```jldoctest
julia> pos, vel, ticks, delta, block = get_slime_bounces(0.0, 1.0, (100, 100), 1e-4);

julia> [pos vel ticks delta block][sortperm(delta, by=abs), :]
42×5 Matrix{Any}:
 -36.3387    (-36.3387, 8.55311)   (48, 99)   2.72403e-6  -36.1875
 -35.3387    (-35.3387, -44.8919)  (98, 49)   2.72403e-6  -35.1875
   4.91125   (4.91125, -9.99654)   (71, 33)  -3.43297e-6    5.0625
   3.59876   (3.59876, -5.95435)   (67, 49)   8.31251e-6    3.75
[...]

julia> sim_bounces(0.0, 1.0, (48, 99))
(2, 1.9600000102072954, 0.9204000199310483)
(3, 2.8804000301383437, 0.8811960291799086)
(4, 3.7615960593182525, 0.842384077962402)
[...]
```
"""
function get_slime_bounces(pos::Float64, vel::Float64, ticks::NTuple{N, Int}, threshold::Float64; entity_type=Pearl, explosion_height=0.061250001192092896, min_pos=pos, max_pos=256., max_ticks=sum(ticks)) where N
  pos, vel = tick_y(entity_type, pos, vel)
  return _get_slime_bounces(BounceIndices(ticks, pos, vel, entity_type), threshold, entity_type, explosion_height, min_pos, max_pos, max_ticks)
end

function _get_slime_bounces(iter::BounceIndices{N}, threshold::Float64, entity_type::Type, explosion_height::Float64, min_pos::Float64, max_pos::Float64, max_ticks::Int) where N
  outpos = Float64[]
  outvel = Float64[]
  outticks = NTuple{N, Int}[]
  outdelta = Float64[]
  outblock = Float64[]
  for BI ∈ iter
    pos = first(BI.P)
    if pos < min_pos || pos > max_pos || sum(BI.I) > max_ticks
      continue
    end
    block, delta = _nearestblockpos(entity_type, pos, explosion_height)
    if abs(delta) < threshold
      push!(outpos, pos)
      push!(outvel, first(BI.V))
      push!(outticks, reverse(BI.I) .- 1)
      push!(outdelta, delta)
      push!(outblock, block)
    end
  end
  sp = sortperm(outdelta, by=abs)
  (outpos[sp], outvel[sp], outticks[sp], outdelta[sp], outblock[sp])
end

"""
    sim_bounces(pos::Float64, indices::Tuple; entitytype::Type)

Simulates a set of bounces listed in `indices`, starting at `pos`;

Useful for reading the output of [`get_slime_bounces`](@ref).

See also [`get_slime_bounces`](@ref).

# Examples
```jldoctest
julia> sim_bounces(0.0, (45, 57))
Bounced: +0.51, long pulse
(0, 0.51, 1.0)
(1, 1.51, 0.9600000102072954)
(2, 2.470000010207295, 0.9204000199310483)
(3, 3.3904000301383435, 0.8811960291799086)
[...]
(45, 11.035809967136265, -0.45525796911456506)
Bounced: +1.0, long pulse
(45, 12.035809967136265, 1.0)
(46, 13.035809967136265, 0.9600000102072954)
[...]
(101, 16.1953340411928, -0.7215951635213262)
(102, 15.473738877671474, -0.7443792180972285)
```
"""
function sim_bounces(pos::Float64, vel::Float64, indices::Tuple; entitytype::Type=Pearl)
  tick = 1
  println((0, pos, vel))
  pos, vel = tick_y(entitytype, pos, vel)
  println((1, pos, vel))
  pos, vel = tick_y(entitytype, pos, vel)
  for (i, bounce) in enumerate(indices)
    for _ in 2:bounce
      println((tick += 1, pos, vel))
      pos, vel = tick_y(entitytype, pos, vel)
    end
    if i != length(indices)
      println("Piston extension on tick $(tick + 1)")
      for d in _slime_bounce(entitytype, pos)
        println((tick += 1, pos += d, 1e0))
        pos, vel = tick_y(entitytype, pos, 1e0)
      end
      println("Piston retraction on tick $(tick + 1)")
    end
  end
  println((tick += 1, pos, vel))
end

function _slime_bounce(::Type{PearlOld}, pos::Float64)
  pos -= floor(pos)
  return 0.75 <= pos ? (0.51,) :
    0.5 <= pos ? (0.,) :
    0.25 <= pos ? (0., 0.51) :
     (0., 0.)
end

function _slime_bounce(::Type{Pearl}, pos::Float64)
  pos -= floor(pos)
  return 0.75 <= pos ? (0.51,) :
    # 1.5 <= pos + 0.9603000092506409 ? (0.,) :
    0.5 <= pos ? (0.,) : # choice between (0) and (0,0.51)
    1.25 <= pos + 0.9603000092506409 ? (0., 0.51) :
     (0., 0.)
end

function _slime_bounce_pos(::Type{PearlOld}, pos::Float64)
  pos -= floor(pos)
  return 0.75 <= pos ? 0.51 :
    0.5 <= pos ? 0. :
    0.25 <= pos ? 1.51 :
    1.
end

function _slime_bounce_pos(::Type{Pearl}, pos::Float64)
  pos -= floor(pos)
  return 0.75 <= pos ? 0.51 :
    0.5 <= pos ? 0. :
    1.25 <= pos + 0.9603000092506409 ? 0.9603000092506409 + 0.51 :
    0.9603000092506409
end

