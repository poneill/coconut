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

Let's also define some type synonyms for clarity and mnemonic ease:

> type BoundFraction = Float
> type Time = Float
> type Force = Float
> type Mass = Float
> type Acceleration = Float

We begin the modeling proper by defining the ligand concentration as a
function of position.  For now, let us assume that that the ligand
concentration is a simple linear function of displacement in the y axis:

> l :: Vector -> Float
> l (Vector x y) = y

We now define the displacement of the cluster and its associated
nodes.  Let's assume for the sake of concreteness that the cluster
possesses 8 receptors spaced with radial symmetry about the cluster
boundary.  

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

> f :: BoundFraction -> Time -> Force
> f bf t = bf * (1 - cos(bf * t))