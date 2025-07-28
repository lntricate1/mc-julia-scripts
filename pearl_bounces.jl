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

function testnew_impl(N::Int, n::Int, M::Int)
  if n == 1
    return quote
      # p1 = p2
      p1 = p2 + _slime_bounce_pos(Pearl, p2)
      v1 = 1.0
      for t1 in 1:$M
        # println(($((Symbol("p", j) for j in 1:N)...),), "; ", ($((Symbol("t", j) for j in 1:N)...),))
        p1, v1 = tick_y(Pearl, p1, v1)
        # println(p1)

        block, delta = _nearestblockpos(Pearl, p1, explosion_height)
        if abs(delta) < threshold
          if block - floor(block) == 0.5
          println(($((Symbol("p", j) for j in 1:N)...),))
          end
          push!(outpos, p1)
          push!(outvel, v1)
          push!(outticks, ($((Symbol("t", j) for j in 1:N)...),))
          push!(outdelta, delta)
          push!(outblock, block)
        end
      end
    end
  else
    tn = Symbol("t", n)
    pn = Symbol("p", n)
    vn = Symbol("v", n)
    pn1 = Symbol("p", n+1)
    vn1 = Symbol("v", n+1)
    return n == N ? quote
      $pn = pos
      $vn = vel
      for $tn in 1:$M
        $(testnew_impl(N, n-1, M))
        $pn, $vn = tick_y(Pearl, $pn, $vn)
      end
    end : quote
      $pn = $pn1 + _slime_bounce_pos(Pearl, $pn1)
      $vn = 1.0
      # $vn = (1-0.03)*0.99f0
      for $tn in 1:$M
        $(testnew_impl(N, n-1, M))
        $pn, $vn = tick_y(Pearl, $pn, $vn)
      end
    end
  end
end

@generated function testnew(pos::Float64, vel::Float64, bounces::Val{N}, ticks::Val{M}, threshold=0.0, explosion_height::Float64=0.061250001192092896) where {N, M}
  # pn = Symbol("p", N+1)
  # vn = Symbol("v", N+1)
  return quote
    outpos = Float64[]
    outvel = Float64[]
    outticks = NTuple{N, Int}[]
    outdelta = Float64[]
    outblock = Float64[]

    # $pn, $vn = tick_y(Pearl, pos, vel)
    $(testnew_impl(N, N, M))
    sp = sortperm(outdelta, by=abs)
    (outpos[sp], outvel[sp], outticks[sp], outdelta[sp], outblock[sp])
    # nothing
  end
end

struct BounceIndices{N,T}
  indices::NTuple{N, Int}
  pos::Float64
  vel::Float64
  entity_type::Type{T}
end

struct BounceIndex{N}
  I::NTuple{N, Int}
  P::NTuple{N, Float64}
  V::NTuple{N, Float64}
end

@inline function __inc(state::Tuple{Int}, pos::Tuple{Float64}, vel::Tuple{Float64}, indices::Tuple{Int}, entity_type::Type{T}) where T
  state_ = state[1] + 1
  pos_, vel_ = tick_y(T, pos[1], vel[1])
  return (state_,), (pos_,), (vel_,)
end

@inline function __inc(state::Tuple{Int, Int, Vararg{Int}}, pos::Tuple{Float64, Float64, Vararg{Float64}}, vel::Tuple{Float64, Float64, Vararg{Float64}}, indices::Tuple{Int, Int, Vararg{Int}}, entity_type::Type{T}) where T
  if state[1] == indices[1]
    state_, pos_, vel_ = __inc(Base.tail(state), Base.tail(pos), Base.tail(vel), Base.tail(indices), T)
    pos__, vel__ = tick_y(T, pos_[1], 1.0)
    return (1, state_...), (pos__ + _slime_bounce_pos(T, pos_[1]), pos_...), (vel__, vel_...)
  end
  pos_, vel_ = tick_y(entity_type, pos[1], vel[1])
  return (state[1] + 1, Base.tail(state)...), (pos_, Base.tail(pos)...), (vel_, Base.tail(vel)...)
end

# @inline function __inc(state::NTuple{N, Int}, pos::NTuple{N, Float64}, vel::NTuple{N, Float64}, indices::NTuple{N, Int}, entity_type::Type) where N
#   I = first(state)
#   ts, tp, tv = Base.tail(state), Base.tail(pos), Base.tail(vel)
#   if I < first(indices)
#     NP, NV = tick_y(entity_type, first(pos), first(vel))
#     return true, (I + 1, ts...), (NP, tp...), (NV, tv...)
#   end
#   first_zero, I, P, V = __inc(ts, tp, tv, Base.tail(indices), entity_type)
#   FP = first(P)
#   NP = first_zero ? FP + _slime_bounce_pos(entity_type, FP) : FP
#   return false, (1, I...), (NP, P...), (1., V...)
# end

@inline function __first(pos::Tuple{Float64}, vel::Float64, entity_type::Type{T}) where T
  p, v = tick_y(T, pos[1], vel)
  return (p,), (v,)
end

@inline function __first(pos::Tuple{Float64, Float64, Vararg{Float64}}, vel::Float64, entity_type::Type{T}) where T
  p1, v1 = __first(Base.tail(pos), vel, T)
  p, v = tick_y(T, p1[1], 1.0)
  return (p + _slime_bounce_pos(T, p), p1...), (v, v1...)
end

@inline function Base.iterate(iter::BounceIndices{N}) where N
  # pos, vel = tick_y(iter.entity_type, iter.pos, iter.vel)
  # BI = BounceIndex(ntuple(i -> 1, Val(N)), ntuple(i -> pos, Val(N)), ntuple(i -> vel, Val(N)))
  BI = BounceIndex(ntuple(i -> 1, Val(N)), __first(ntuple(i -> iter.pos, Val(N)), iter.vel, iter.entity_type)...)
  return BI, BI
end

@inline function Base.iterate(iter::BounceIndices{N,T}, state::BounceIndex{N}) where {N,T}
  state.I == iter.indices && return nothing
  # _, I, P, V = __inc(state.I, state.P, state.V, iter.indices, iter.entity_type)
  I, P, V = __inc(state.I, state.P, state.V, iter.indices, iter.entity_type)
  next = BounceIndex(I, P, V)
  # I = __inc(state.I, iter.indices, iter.entity_type)
  # next = BounceIndex(I, state.P, state.V)
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
  # pos, vel = tick_y(entity_type, pos, vel)
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
function sim_bounces(pos::Float64, vel::Float64, indices::Tuple; entity_type::Type=Pearl)
  return sim_bounces_old(pos, vel, indices .+ 1; entitytype=entity_type)
end

function sim_bounces_old(pos::Float64, vel::Float64, indices::Tuple; entitytype::Type=Pearl)
  tick = 0#1
  println((0, pos, vel))
  pos, vel = tick_y(entitytype, pos, vel)
  # println((1, pos, vel))
  # pos, vel = tick_y(entitytype, pos, vel)
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
  println("_slime_bounce(Pearl, $pos)")
  return 0.75 <= pos ? (0.51,) :
    # 1.5 <= pos + 0.97*0.99f0 ? (0.,) :
    0.5 <= pos ? (0.,) : # choice between (0) and (0,0.51)
    1.25 <= pos + 0.97*0.99f0 ? (0., 0.51) :
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
    1.25 <= pos + 0.97*0.99f0 ? 0.97*0.99f0 + 0.51 :
    0.97*0.99f0
end
