This file documents the first prototype of a model for border cell
migration in Drosophila.  We elect to model the border cell cluster as
a rigid homogenous spherical mass in a ligand gradient.  The cluster
possesses ligand receptor nodes on its boundary which exert forces
proportional to the bound ligand fraction and normal to the surface.
Additionally, we assume that the motive force is driven by a forcing
function that models the saturation/desensitization dynamics of the
receptor.

Let us begin by sketching a vector type which will allow us to
simplify expressions later on.

> data Vector = Vector Float Float deriving (Eq, Show)

We define addition,

> (.+) :: Vector -> Vector -> Vector
> (.+) (Vector a b) (Vector c d) = Vector (a + c) (b + d)

subtraction,

> (.-) :: Vector -> Vector -> Vector
> (.-) (Vector a b) (Vector c d) = Vector (a - c) (b - d)

and scalar multiplication.

> (.*) :: Float -> Vector -> Vector
> (.*) c (Vector a b) = Vector (c * a) (c * b)

Let's also define a Node type for convenience.  A node is just a
position:

> type Index = Int
> data Node = Node Vector
> type Nodes = [Node]

We will often be referring to the center of the cluster, and may wish
to designate it by a type synonym:  

> type Center = Vector

Let's also define some type synonyms for clarity and mnemonic ease:

> type BoundFraction = Float
> type Time = Float
> type Force = Vector
> type Mass = Float
> type Acceleration = Vector

We begin the modeling proper by defining the ligand concentration as a
function of position.  For now, let us assume that that the ligand
concentration is a simple linear function of displacement in the y axis:

> l :: Vector -> Float
> l (Vector x y) = y

We now define the displacement of the cluster and its associated
nodes.  Let's assume for the sake of concreteness that the cluster
possesses 8 receptors spaced with radial symmetry about the cluster
boundary.  

> numNodes = 8
> nodeList = [0..numNodes - 1]

The function theta describes the normal vector for the ith node:

> theta :: Index -> Vector
> theta i = Vector (cos (2 * pi * i' / n)) (sin (2 * pi * i' / n))
>   where i' = fromIntegral i 
>         n = fromIntegral numNodes          

the position of the ith node is just the position of the center,
displaced by the length of the radius in the direction of the ith normal vector:

> ithNode :: Center -> Index -> Vector 
> ithNode c i = c .+ (radius .* (theta i))

where

> radius = 1

For now, let the total receptor capacity of the ith node be constant: 

> rt :: Int -> Float
> rt i = 10

Next, let's define the fraction of ligand bound to a receptor:

> ri :: Vector -> Index -> BoundFraction
> ri v i = (lv * (rt i))/(lv + k) 
>   where lv = l v 

The above expression requires a definition for the
equilibrium constant K, so let

> k = 1

Next we define a forcing function f which describes the force exerted
by the ith node as a function of time:

> timestep = 0.001
> f :: BoundFraction -> Time -> Float
> f bf t = bf * (1 - cos(bf * t * timestep))

> ithForce :: Center -> Index -> Time -> Vector
> ithForce c i t = (f (ri v i) t) .* theta i 
>   where v = ithNode c i
  
Now we are in a position to sum the forces acting on each node in
order to obtain an expression for the force acting on the center of
the cluster:

> forceOnCenter :: Center -> Time -> Force
> forceOnCenter c t = foldl1 (.+) [forceFrom i | i <- nodeList]
>   where forceFrom = \i -> ithForce c i t


Finally we can write the expression for the acceleration on the center
at time t:

> pos'' :: Center -> Time -> Acceleration
> pos'' c t = (1/m) .* (forceOnCenter c t)

where m is the mass of the cluster:

> m = 1

Lastly let us write a simple step function that updates the
displacement of the center of the cluster over a small timestep.

We should read the following signature as a function that accepts a
Center vector (corresponding to the initial condition), a time (at
which we wish to know the position of the center) and a timestep, and
returns a center vector describing the position of the center at the
desired time.

> posRef :: Center -> Time -> Float -> Center
> posRef c t dt 
>   | t == 0 = c
>   | t == dt = c
>   | otherwise  =  posTerm .+ accTerm
>   where s   = \t -> posRef c t dt
>         s'' = pos'' c
>         posTerm = (2 .* s (t - dt)) .- (s (t - 2 * dt))
>         accTerm = ((dt**2) .* (s'' (t - 2 * dt)))

Or we can think of pos as a function which accepts the positions at
the last two time steps and returns the position at the current
timestep.  Note that the timestep is implicitly given by the indices
of the list of vectors.  Although we end up discarding almost all of
the data, the only alternative would be to track the index explicitly
with another parameter, say by zipping the displacement history with
the list [0..]

> pos :: [Center] -> Center 
> pos cs = posTerm .+ accTerm
>   where c = head cs
>         s'' = pos'' c
>         n = length cs
>         [oneBack, twoBack] = drop (n - 2) cs
>         twoAgo = fromIntegral n - 2
>         posTerm = (2 .* oneBack) .- twoBack
>         accTerm = (1 .* (s'' twoAgo))
