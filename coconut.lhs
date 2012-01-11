\documentclass{article}
\usepackage{fancyvrb}
\DefineVerbatimEnvironment{code}{Verbatim}{fontsize=\small}
\DefineVerbatimEnvironment{example}{Verbatim}{fontsize=\small}
\newcommand{\ignore}[1]{}
\usepackage{verbatim}
\title{A model for border cell migration in \textit{Drosophila}}
\author{Patrick O'Neil}
\begin{document}
\maketitle{}
This file documents the first prototype of a model for border cell
migration in \textit{Drosophila}.  It is both an executable \LaTeX
file which describes the model, and an executable Haskell file that
implements the model.

\section{Background}
During oogenesis in \textit{Drosophila}, an egg chamber develops from
a set of 16 germ cells surrounded by approximately 1,000 follicle
cells.  One of the germ cells differentiates into the oocyte and grows
to roughly half the total volume of the egg chamber. The other germ
cells are fated to become \textit{nurse cells}, polyploid cells whose
ultimate purpose will be to contribute nutrients to the maturing
oocyte.  During stage 9 of oocyte development, two \textit{polar
  cells} on the anterior tip of the oocyte induce differentiation of
neighboring cells into \textit{border cells}, forming a polar
cell-border cell complex referred to as the \textit{border cell
  cluster}.  The cluster then delaminates from the epithelium and
traverses the egg chamber in order to reach the oocyte.  The motility
of the border cell complex is governed by three chemoattractant
molecules expressed by the oocyte, which bind to tyrosine kinase
receptors in the border cells and induce actin polymerization,
generating protrusions which gain purchase against the tightly packed
nurse cells through adherens junctions.  Despite the importance of
border cell migration in \textit{Drosophila} as a general model for
cell migration and invasion, little is understood about the specific
mechanisms of cell motility.  In particular, several questions have
arisen from observations that the cluster appears to rotate or spin as
it negotiates the nurse cell gauntlet:
\begin{enumerate}
\item What can account for the rotation of the cluster?
\item Can a model which assumes independence among the border cells
  explain the observed behavior, or must their efforts be coordinated
  through a signaling mechanism?
\end{enumerate}
 We elect to
model the border cell cluster as a rigid homogenous spherical mass in
a ligand gradient.  The cluster possesses ligand receptor nodes on its
boundary which exert forces proportional to the bound ligand fraction
and normal to the surface.  Additionally, we assume that the motive
force is driven by a forcing function that models the
saturation/desensitization dynamics of the receptor.

\section{Data types and synonyms}
In this section we declare some data types and type synonyms that will
aid us later on.  In general, the significance of the code in this
section should be clear from context later on.  However, the reader
may wish to note that operators preceded by a \texttt{.}
(\textit{e.g.} \texttt{(.+)}) just represent the relevant vector
operations.  This is a quirk of Haskell's type system which the reader
is urged to ignore.

Let us begin by sketching a vector type which will allow us to
simplify expressions later on.

\begin{code}
import Debug.Trace
\end{code}

\begin{code}
data Vector = Vector Float Float deriving (Eq, Show)
\end{code}

We define addition,

\begin{code}
(.+) :: Vector -> Vector -> Vector 
(.+) (Vector a b) (Vector c d) = Vector (a + c) (b + d)
\end{code}

subtraction,

\begin{code}
(.-) :: Vector -> Vector -> Vector 
(.-) (Vector a b) (Vector c d) = Vector (a - c) (b - d)
\end{code}

and (left) scalar multiplication.

\begin{code}
(.*) :: Float -> Vector -> Vector 
(.*) c (Vector a b) = Vector (c * a) (c * b)
\end{code}

\begin{code}
norm :: Vector -> Float
norm (Vector x y) = sqrt (x ** 2 + y ** 2)
\end{code}
Let's also define a Node type for convenience.  A node is just a
position:

\begin{code}
type Index = Int 
data Node = Node Vector 
type Nodes = [Node]
\end{code}

We will often be referring to the center of the cluster, and may wish
to designate it by a type synonym:

\begin{code}
type Center = Vector
\end{code}

Let's also define some type synonyms for clarity and mnemonic ease:

\begin{code}
type BoundFraction = Float 
type Time = Float 
type Force = Vector
type Mass = Float 
type Acceleration = Vector
\end{code}

\section{Model Definition}

\subsection{Ligand concentration}
We begin the modeling proper by defining the ligand concentration as a
function of position.  For now, let us assume that that the ligand
concentration is a simple function of distance from the point (10,10).
This clearly biologically implausible, but will be tolerated for now
in the interest in getting the model up and running.

\begin{code}
emitterLocation = Vector 10 10
l :: Vector -> Float
l (Vector x y) = 1 / (d ** 2)
                 where Vector x' y' = emitterLocation
                       d = sqrt ((x - x')**2 + (y - y')**2)
\end{code}

\subsection{Cluster configuration}
We now define the displacement of the cluster and its associated
nodes.  Let's assume for the sake of concreteness that the cluster
possesses 8 receptors spaced with radial symmetry about the cluster
boundary.  

\begin{code}
numNodes = 8
nodeList = [0..numNodes - 1]
\end{code}

The function theta gives the normal unit vector for the ith node:

\begin{code}
theta :: Index -> Vector
theta i = Vector (cos (2 * pi * i' / n)) (sin (2 * pi * i' / n))
  where i' = fromIntegral i 
        n = fromIntegral numNodes          
\end{code}

The position of the ith node is just the position of the center,
displaced by the length of the radius in the direction of the ith normal vector:

\begin{code}
ithNode :: Center -> Index -> Vector 
ithNode c i = c .+ (radius .* (theta i))
\end{code}

where

\begin{code}
radius = 1
\end{code}

\subsection{Ligand-receptor binding}
For now, let the total receptor capacity of the ith node be constant: 

\begin{code}
rt :: Int -> Float
rt i = 1
\end{code}

Next, let's define the fraction of ligand bound to a receptor:

\begin{code}
ri :: Vector -> Index -> BoundFraction
ri v i = (lv * (rt i))/(lv + k) 
  where lv = l v 
\end{code}

The above expression requires a definition for the equilibrium
constant K, so we arbitrarily let

\begin{code}
k = 1
\end{code}

Next we define a forcing function f which describes the force exerted
by the ith node as a function of time:

\begin{code}
timestep = 1
f :: BoundFraction -> Time -> Float
f bf t = bf * (1 - cos(bf * t * timestep))
\end{code}

\begin{code}
ithForce :: Center -> Index -> Time -> Vector
ithForce c i t    = magnitude .* direction
  where v         = ithNode c i
        magnitude = (f (ri v i) t)
        direction = theta i 
\end{code}

\subsection{Force balancing} 
Now we are in a position to sum the forces acting on each node in
order to obtain an expression for the force acting on the center of
the cluster:

\begin{code}
forceOnCenter :: Center -> Time -> Force
forceOnCenter c t = foldl1 (.+) [forceFrom i | i <- nodeList]
  where forceFrom = \i -> ithForce c i t
\end{code}


Finally we can write the expression for the acceleration on the center
at time \textit{t}:

\begin{code}
pos'' :: Center -> Time -> Acceleration
pos'' c t = (1/m) .* (forceOnCenter c t)
\end{code}

where m is the mass of the cluster:

\begin{code}
m = 100
\end{code}

\subsection{Iteration}
Lastly let us write a simple step function that updates the
displacement of the center of the cluster over a small timestep.
Consider that we must describe the position of the cluster $s(t)$ as a
function of time, given the acceleration function $s''(t)$.  We assume
for now that the displacement of the cluster is determined only by the
acceleration profile and the initial position.  We assume no damping
due to friction or viscosity.  

To motivate our approach, recall that we may write:

$$s(t + \Delta t) \approx s(t) + s'(t)\Delta t$$

for $\Delta t$ sufficiently small.  This, alas, does not help since we
do not know $s'(t)$ in general.  However we may appeal to the
definition of the difference quotient and write:

$$s'(t) \approx \frac{s(t + \Delta t) - s(t)}{\Delta t}.$$

Substituting and simplifying, we arrive at:

$$s(t + \Delta t) \approx 2s(t) - s(t - \Delta t) + s''(t - \Delta t)\Delta^2t$$

which gives us the $t + \Delta t$ position after the time step as a
function of the $t$th and $t - \Delta t$th positions and the $t +
\Delta t$th acceleration value.  

Let's summarize the model so far.  We have described the acceleration
on the cluster as a function of position and time.  Formally, we may write:

$$a(\mathbf{x},t) = \frac{1}{m}\displaystyle\sum{i=1}^n\mathbf{g}(\mathbf{u}_i,t)$$

where $m$ is the cluster mass and $\mathbf{g}(\mathbf{u}_i,t)$ is the force due
to the $i$th node at time $t$.  In turn, $g$ is defined as:

$$\mathbf{g}(r_i(\mathbf{u}),t) = f(\mathbf{u},t)\mathbf{n}$$

where $r_i$ gives the bound fraction of ligand at the point
$\mathbf{u}$, $\mathbf{n}$ is the normal vector to the cluster at
$\mathbf{u}$, and $f$ is a forcing function meant to model the
periodic protrusion and relaxation of the polymerized actin.  If we
assume that the reaction is at equilibrium then, after some
manipulation of the mass action law, we can express the bound fraction
at position $\mathbf{u}$ as

$$r_i(\mathbf{u}) = \frac{L(\mathbf{u})R_{it}}{L(\mathbf{u}) + K}$$

  where $R_{it}$ is the receptor total at the node, and $K$ is the
  equilibrium constant for the ligand binding reaction.  

  The unit normal vector $\mathbf{n}$ is just:

  $$\mathbf{n}(\mathbf{u}) = \frac{\mathbf{u} - \mathbf{c}}{R}$$
  
  where $\mathbf{c}$ is the position of the center of the cluster and
  $R$ is its radius.
\begin{code}
s'' t = (-10)
dt = 0.01
s t | trace (show t) False = undefined
s t 
  | t <= 0 = 100
  | otherwise = 2 * (s (t - dt)) - (s (t - 2 * dt)) + (s'' (t - dt)) * dt ** 2
\end{code}

\begin{code}
q'' t = (-10)
q t1 t2 = 2 * t1 - t2 + (q'' t1) * dt**2
\end{code}

We should read the following signature as a function that accepts a
Center vector (corresponding to the initial condition), a time (at
which we wish to know the position of the center) and a timestep, and
returns a center vector describing the position of the center at the
desired time.

\begin{code}
posRef :: Center -> Time -> Float -> Center
posRef c t dt 
  | t == 0 = c
  | t == dt = c
  | otherwise  =  posTerm .+ accTerm
  where s   = \t -> posRef c t dt
        s'' = pos'' c
        posTerm = (2 .* s (t - dt)) .- (s (t - 2 * dt))
        accTerm = ((dt**2) .* (s'' (t - 2 * dt)))
\end{code}

Or we can think of pos as a function which accepts the positions at
the last two time steps and returns the position at the current
timestep.  Note that the timestep is implicitly given by the indices
of the list of vectors.  Although we end up discarding almost all of
the data, the only alternative would be to track the index explicitly
with another parameter, say by zipping the displacement history with
the list [0..]

\begin{code}
pos :: [Center] -> [Center]
pos cs = cs ++ [posTerm .+ accTerm]
  where c = head cs
        s'' = pos'' c
        n = length cs
        [twoBack, oneBack] = drop (n - 2) cs
        twoAgo = fromIntegral n - 2
        posTerm = (2 .* oneBack) .- twoBack
        accTerm = (1 .* (s'' twoAgo))
\end{code}


step :: Center -> [Center]
step c = iterate pos' [c, c]
  where pos' cs = cs ++ [pos cs]

\end{document}