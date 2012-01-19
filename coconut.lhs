\documentclass{article}
\usepackage{fancyvrb}
\usepackage{graphicx}
\usepackage{amsmath}
\DefineVerbatimEnvironment{code}{Verbatim}{fontsize=\small}
\DefineVerbatimEnvironment{example}{Verbatim}{fontsize=\small}
\newcommand{\ignore}[1]{}
\newcommand{\deltat}{\Delta t}
\newcommand{\deltatwo}{\frac{\Delta t}{2}}
\newcommand{\deltafour}{\frac{\Delta t}{4}}
\usepackage{verbatim}
\title{Sketch of a model for border cell migration in \textit{Drosophila}}
\begin{document}
\maketitle{}
This file documents the first prototype of a model for border cell
migration in \textit{Drosophila}.  It is both an executable \LaTeX 
file which describes the model, and an executable Haskell file that
implements the model.

\section{Background}
In lieu of a reference section, let me mention that everything in this section
(except the mistakes) is cribbed from various papers by MSG and
Montell lab.  

During oogenesis in \textit{Drosophila}, an egg chamber develops from
a set of 16 germ cells surrounded by approximately 1,000 follicle
cells.  One of the germ cells differentiates into the oocyte and grows
to roughly half the total volume of the egg chamber. The other germ
cells are fated to become \textit{nurse cells}, polyploid cells whose
ultimate purpose will be to contribute ngantletutrients to the
maturing oocyte.  During stage 9 of oocyte development, two
\textit{polar cells} on the anterior tip of the oocyte induce
differentiation of neighboring cells into \textit{border cells},
forming a polar cell-border cell complex referred to as the
\textit{border cell cluster}.  The cluster then delaminates from the
epithelium and traverses the egg chamber in order to reach the oocyte.
The motility of the border cell complex is governed by three
chemoattractant molecules expressed by the oocyte which bind to
tyrosine kinase receptors in the border cells and induce actin
polymerization, generating protrusions which gain purchase against the
tightly packed nurse cells via adherens junctions.  Despite the
importance of border cell migration in \textit{Drosophila} as a
general model for cell migration and invasion events, little is
understood about the specific mechanisms of cell motility.  In
particular, several questions have arisen from observations that the
cluster appears to rotate or spin as it negotiates the nurse cell
``gantlet''.
\begin{enumerate}
\item What can account for the rotation of the cluster?
\item Can a model which assumes independence among the border cells
  explain the observed behavior, or must their efforts be coordinated
  through a signaling mechanism?
\end{enumerate}
We first elect to model the border cell cluster as a rigid homogenous
spherical mass in a ligand gradient.  The cluster possesses ligand
receptor nodes on its boundary which exert forces proportional to the
bound ligand fraction and normal to the surface.  Additionally, we
assume that the motive force is driven by a forcing function that
models the saturation/desensitization dynamics of the receptor.

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
data Vector = Vector { getX :: Float 
                     , getY :: Float
                     } deriving (Eq, Show)
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

We might as well have an inner product as well:
\begin{code}
(.*.) :: Vector -> Vector -> Float
(.*.) (Vector x y) (Vector x' y') = x * x' + y * y'
\end{code}
\begin{code}
norm :: Vector -> Float
norm v = sqrt $ v .*. v
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
This is biologically doubtful, but will be tolerated for now
in the interest in getting the model up and running.

\begin{code}
emitterLocation = Vector 20 20
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
ithNode c i = c .+ (radius .* theta i)
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
ri v i = (lv * rt i)/(lv + k) 
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
--f bf t = bf * (1 - cos(bf * t * timestep)) debugging
f bf t = bf
\end{code}

\begin{code}
ithForce :: Center -> Index -> Time -> Vector
ithForce c i t    = magnitude .* direction
  where v         = ithNode c i
        magnitude = f (ri v i) t
        direction = theta i 
\end{code}

\subsection{Force balancing} 
Now we are in a position to sum the forces acting on each node in
order to obtain an expression for the force acting on the center of
the cluster:

\begin{code}
forceOnCenter :: Center -> Time -> Force
forceOnCenter c t = foldl1 (.+) [forceFrom i | i <- nodeList]
  where forceFrom i = ithForce c i t
\end{code}


Finally we can write the expression for the acceleration on the center
at time \textit{t}:

\begin{code}
pos'' :: Center -> Time -> Acceleration
pos'' c t = (1/m) .* forceOnCenter c t
\end{code}

where m is the mass of the cluster:

\begin{code}
m = 1
\end{code}


\subsection{Numerical Simulation}
Lastly let us write a  step function that updates the
displacement of the center of the cluster over a small timestep.
Consider that we must describe the position of the cluster $s(t)$ as a
function of time, given the acceleration function $s''(t)$.  We assume
for now that the displacement of the cluster is determined only by the
acceleration profile and the initial position.  We assume no damping
due to friction or viscosity.  

To motivate our approach, recall that we may write down the following
difference equation due to Euler:

$$s(t + \Delta t) \approx s(t) + s'(t)\Delta t$$

for $\Delta t$ sufficiently small.  One may observe that in the Euler
method, the $t + \deltat$th position update is informed by the
derivative at time $t$, which in general is a poorer guide to the
behavior of the function over the interval $[t, t + \deltat]$ than the
derivative at the midpoint of the interval, $t + \deltatwo$.  That is,
we would like to write instead: 

$$s(t + \Delta t) \approx s(t) + s'(t + \deltatwo)\Delta t$$.

Here, we may fear the spectre of infinite regress, since the
computation of $s(t + \deltat)$ requires $s'(t + \deltatwo)$, which in
turn would seem to require division of the timestep by increasingly
large powers of two.  Instead, we approximate $s'(t + \deltatwo)$ via
Taylor expansion.  Recall that we may write any (sufficiently
well-behaved) function $f$ as a power series expanded about an
arbitrary point $a$\footnote{Provided a positive radius of convergence and the containment of $x$ by that radius, but these details will be happily neglected.}:

$$f(x) = f(a) + f'(a)(x - a) + \frac{f''(a)(x - a)^2}{2} + \ldots + \frac{f^{(n)}(a)(x - a)^n}{n!}$$

Taylor-expanding $s'$ about $t$ out to one term, we have:

$$s'(t + \deltatwo) = s'(t) + \deltatwo s''(t)$$
and substitution back into the original difference equation yields:
$$s(t + \Delta t) \approx s(t) + \deltat [s'(t) + \deltatwo s''(t)].$$

Finally, we can replace $s'(t)$ by the definition of the difference quotient and write:

$$s(t + \Delta t) \approx s(t) + \deltat [[\frac{s(t) - s(t - \deltat)}{\deltat} + \deltatwo s''(t)]$$

which simplifies to:

$$s(t + \Delta t) \approx 2s(t)  - s(t - \deltat) + \frac{\Delta^2 t}{2} s''(t)]$$


giving us the $t + \Delta t$ position after the time step as a
function of the $t$th and $t - \Delta t$th positions and the $t$th
acceleration value.

\subsection{Summary of Model Definition}
Let's summarize the model so far.  We have described the acceleration
on the cluster as a function of position and time.  Formally, we may write:

$$a(\mathbf{x},t) = \frac{1}{m}\displaystyle\sum_{i=1}^n\mathbf{g}(\mathbf{u}_i,t)$$

where $m$ is the cluster mass and $\mathbf{g}(\mathbf{u}_i,t)$ is the force due
to the $i$th node at time $t$.\footnote{(The function name $f$ having already been taken.)}  In turn, $g$ is defined as:

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

Below is a straightforward implementation of (Eq \ref{eq:step}):

\begin{example}
pos :: Center -> Time -> Float -> Center
pos c t dt 
  | t == 0 = c
  | t == dt = c
  | otherwise  =  posTerm .+ accTerm
  where s   = \t -> posRef c t dt
        s'' = pos'' c
        posTerm = (2 .* s (t - dt)) .- (s (t - 2 * dt))
        accTerm = ((dt**2/2) .* (s'' (t - 1 * dt))) --NB computing t from t - 1
\end{example}

This implementation, however, suffers exponential slowdown (since $s(t
- 2 \Delta t)$ will be recomputed from scratch during the computation
of $s(t - \Delta t)$).  More efficiently, we can think of \texttt{pos}
as a function which accepts the entire position history as a list and
returns a new history with the current position appended.  Note that
the timestep is implicitly given by the indices of the list of
vectors.  

\begin{code}
pos :: [Center] -> [Center]
pos cs = cs ++ [posTerm .+ accTerm]
  where c = head cs
--        s'' = pos'' c
        n = length cs
        [twoBack, oneBack] = drop (n - 2) cs
        oneAgo = fromIntegral n - 1
        posTerm = (2 .* oneBack) .- twoBack
        accTerm = (1/2) .* pos'' oneBack oneAgo
\end{code}

While there are further gains in optimization to be had here, they
need not detain us at the moment.

\section{Model performance}

First let's address the first question posed in the introduction and
ask whether the current model is capable of exhibiting rotation.
Recall from the last section that the acceleration on the cluster is
the sum of the accelerations due to each node, each of which are
normal to the surface of the cluster.  If we resolve each acceleration
vector into its normal and tangent components, the tangent component
in every case will be the zero vector, hence no rotation can occur.
Finally, note that this result is independent of the placement of the
nodes on the cluster boundary.

\begin{code}
iterateN :: Int -> (a -> a) -> a -> a 
iterateN 0 f x = x
iterateN n f x = iterateN (n - 1) f (f x)


\end{code}

\begin{code}
c = Vector 0 0
history = iterateN 10000 pos [c, c]
\end{code}

Since the entire model is symmetric with respect to the axis $y = x$,
we may analyze the kinematics with respect to that axis.  Let us
consider the first 5000 iterations and recover the displacements:

\begin{code}
ps =  map norm history
\end{code}

Next we recover an approximation to the velocity profile:

\begin{code}
vs = zipWith (-) (tail ps)  ps
\end{code}

and finally, we can compute the acceleration in two different ways in
order to perform a sanity check:
\begin{code}
as1 = zipWith (-) vs (tail vs) 
as2 = map ( (/ m) . proj . uncurry forceOnCenter) $ zip history [0..]
  where proj = (.*.) (Vector 1 1)
\end{code}

Figure \ref{fig:1} depicts the result of a sample run, in which the
cluster begins at the origin and migrates toward the source of the
chemoattractant at (20,20).  Position, velocity and acceleration are
depicted, with position scaled down by a factor of 10:

\begin{figure}[ht]
  \centering
  \includegraphics[scale=.65]{firstrun.png}
  \caption{Sample run}
  \label{fig:1}
\end{figure}

Several features warrant mention here.  Most strikingly, model does
not reach a steady state in which the cluster is localized to the
chemoattractant source, but instead ``slingshots" the cluster out to
infinity.  This undesirable behavior is no doubt due to the definition
of the ligand concentration function and its singularity at (10,10),
compounded by the numerical error of the simulation method.
Variation of the free parameters (\textit{viz} mass, radius, timestep) do
not seem to change this story significantly.

\section{Future Directions}
The first model seems to leave several features to be desired.  Most
notably, the impossibility of cluster rotation comes as a serious
drawback.  Other considerations can be divided roughly into three
areas of concern: the physical configuration of the cluster, the
nature of the motive force, and the qualities of the surrounding
medium.  

Regarding the first category, one might wish to revise the assumption
of infinite rigidity: movies of the migration event seem to reveal
significant deformation of the cluster as it squeezes through the gaps
between neighboring nurse cells.  One might instead model the cluster
as an elastic body with constant volume or surface area, perhaps
describing the cluster boundary locally as a flexible beam.  To an
even coarser approximation, the boundary might be modeled as a
``necklace" of receptors on springs.  (Conceptually, they might as
well be rigid linkages, but my intuition is that allowing the linkages
to ``give" a bit by introducing a moderate spring constant would alleviate some of the potential for numerical error.  This is just an unfounded hunch.)

Into the second category of concerns fall questions about the
assumption that the protrusions are normal to the cluster boundary.  I
haven't seen much evidence to the contrary in the movies, but it would
be possible to get at least some rotation out of such a model.
Perhaps more to the point is the omission of the adhesion dynamics: by
treating the motive force as a normal vector, we might be missing out
on torques arising from the asymmetric adhesion of protrusions to
nurse cell membranes.

Lastly, the prevailing assumptions about the intercellular medium are
sure candidates for revision in subsequent rounds of model-building.
Without a viscosity/friction term, the medium approximates an ice
rink, and the cluster velocity is unbounded.  In fact, the egg chamber
is probably a better example of ``life at low Reynolds number", in
which inertial forces are dominated by viscous forces.  Perhaps adding
a tuning parameter in the difference equation governing the relative
significance of the acceleration and velocity terms would be the
cleanest way to begin to address that.

\section{Damping}
To address the concern posed in the preceding paragraph concerning the
neglect of viscosity, we modify the equations of motion given in the
last section.  Formally, we have:

\begin{equation}
  \label{eq:dampedAnalytic}
  s''(t) = f(t,s(t)) - k s'(t)
\end{equation}
which captures the intuition that the cluster should experience a
frictional force proportional to its velocity in the direction
opposite its travel.  We may represent this numerically as follows:
\begin{align*}
s(t + \Delta t) \approx& s(t) + s'(t + \deltatwo)\Delta t\\
s'(t + \deltat) \approx& s'(t) + s\prime\prime(t + \deltatwo) \deltat
\label{eq:damped}
\end{align*}
where

$$\prime\prime(t) = f(t,s(t)) - k s'(t)$$

and $k$ has units 1/time.
Let us now begin to cash this out:
\begin{align*}
  s'(t + \deltatwo) \approx& s'(t) + s\prime\prime(t + \deltafour) \deltatwo\\
  \approx& \frac{s(t) - s(t - \deltat)}{\deltat} + [f(t+\deltafour,s(t+\deltafour)) - k s'(t + \deltafour)]\deltatwo\\
\approx& \frac{s(t) - s(t - \deltat)}{\deltat} + [f(t+\deltafour,s(t)+s'(t)\deltafour)) - k(s'(t) + s\prime\prime(t)\deltafour)]\deltatwo\\
\approx& \frac{s(t) - s(t - \deltat)}{\deltat} + [f(t+\deltafour,s(t)+s'(t)\deltafour)) - k(\frac{s(t) - s(t - \deltat)}{\deltat} + s\prime\prime(t)\deltafour)]\deltatwo\\
\approx& \frac{s(t) - s(t - \deltat)}{\deltat} + [f(t+\deltafour,s(t)+\frac{s(t) - s(t - \deltat)}{\deltat}\deltafour)) - k(\frac{s(t) - s(t - \deltat)}{\deltat} + s\prime\prime(t)\deltafour)]\deltatwo\\
\end{align*}
and finally we may write:
$$s(t + \deltat) \approx s(t) + [\frac{s(t) - s(t - \deltat)}{\deltat} + [f(t+\deltafour,s(t)+s'(t)\deltafour)) - k(\frac{s(t) - s(t - \deltat)}{\deltat} + s\prime\prime(t)\deltafour)]\deltatwo]\deltat$$


\end{document}