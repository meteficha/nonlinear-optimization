---------------------------------------------------------------------------
-- | Module    : Math.Statistics.Dirichlet.Mixture
-- Copyright   : (c) 2009 Felipe Lessa
-- License     : GPL
--
-- Maintainer  : felipe.lessa@gmail.com
-- Stability   : experimental
-- Portability : portable
--
-- This module implements the algorithms described by Hager and
-- Zhang [1].  We use bindings to CG_DESCENT library by the same
-- authors, version 3.0 from 18/05/2008 [2].  The library code is
-- also licensed under the terms of the GPL.
--
-- * [1] Hager, W. W. and Zhang, H.  /A new conjugate gradient/
--   /method with guaranteed descent and an efficient line/
--   /search./ Society of Industrial and Applied Mathematics
--   Journal on Optimization, 16 (2005), 170-192.
--
-- * [2] http://www.math.ufl.edu/~hager/papers/CG/CG_DESCENT-C-3.0.tar.gz
--
--------------------------------------------------------------------------


module Math.Optimization.Algorithms.HagerZhang05
    (-- * Options
     Parameters(..)
    ,Verbose(..)
    ,LineSearch(..)
    ,StopRules(..)
    ,EstimateError(..)
     -- * Technical parameters
    ,TechParameters(..)
    ) where

import Foreign
import Foreign.C
#include "cg_user.h"

-- | Parameters given to the optimizer.
data Parameters = Parameters {
    printFinal :: Bool
    -- ^ Print final statistics to @stdout@.  Defaults to @True@.

    ,printParams :: Bool
    -- ^ Print parameters to @stdout@ before starting.  Defaults to @False@

    ,verbose :: Verbose
    -- ^ How verbose we should be while computing.  Everything is
    -- printed to @stdout@. Defaults to 'Quiet'.

    ,lineSearch :: LineSearch
    -- ^ What kind of line search should be used.  Defaults to
    -- @AutoSwitch 1e-3@.

    ,qdecay :: CDouble
    -- ^ Factor in @[0, 1]@ used to compute average cost
    -- magnitude @C_k@ as follows:
    --
    -- > Q_k = 1 + (qdecay)Q_{k-1},   Q_0 = 0
    -- > C_k = C_{k-1} + (|f_k| - C_{k-1})/Q_k
    --
    -- Defaults to @0.7@.

    ,stopRules :: StopRules
    -- ^ Stop rules that define when the iterations should end.
    -- Defaults to @DefaultStopRule 0@.

    ,estimateError :: EstimateError
    -- ^ How to calculate the estimated error in the function
    -- value.  Defaults to @RelativeEpsilon 1e-6@.

    ,quadraticStep :: Maybe CDouble
    -- ^ When to attempt quadratic interpolation in line search.
    -- If @Nothing@ then never try a quadratic interpolation
    -- step.  If @Just cutoff@, then attemp quadratic
    -- interpolation in line search when @|f_{k+1} - f_k| / f_k
    -- <= cutoff@.  Defaults to @Just 1e-12@.

    ,debugTol :: Maybe CDouble
    -- ^ If @Just tol@, then always check that @f_{k+1} - f_k <=
    -- tol * C_k@. Otherwise, if @Nothing@ then no checking of
    -- function values is done.  Defaults to @Nothing@.

    ,initialStep :: Maybe CDouble
    -- ^ If @Just step@, then use @step@ as the initial step of
    -- the line search.  Otherwise, if @Nothing@ then the initial
    -- step is programatically calculated.  Defaults to
    -- @Nothing@.

    ,maxItersFac :: CInt
    -- ^ Defines the maximum number of iterations.  The process
    -- is aborted when @maxItersFac * n@ iterations are done, where
    -- @n@ is the number of dimensions.  Defaults to @maxBound@.

    ,nexpand :: CInt
    -- ^ Maximum number of times the bracketing interval grows or
    -- shrinks in the line search.  Defaults to @50@.

    ,nsecant :: CInt
    -- ^ Maximum number of secant iterations in line search.
    -- Defaults to @50@.

    ,restartFac :: CInt
    -- ^ Restart the conjugate gradient method after @restartFac
    -- * n@ iterations. Defaults to @1@.

    ,funcEpsilon :: CDouble
    -- ^ Stop when @-alpha * dphi0@, the estimated change in
    -- function value, is less than @funcEpsilon * |f|@.
    -- Defaults to @0@.

    ,nanRho :: CDouble
    -- ^ After encountering @NaN@ while calculating the step
    -- length, growth factor when searching for a bracketing
    -- interval.  Defaults to @1.3@.

    ,techParameters :: TechParameters
    -- ^ Technical parameters which you probably should not
    -- touch.
    } deriving (Eq, Ord, Show, Read)

-- | Technical parameters which you probably should not touch.
-- You should read the papers of @CG_DESCENT@ to understand how
-- you can tune these parameters.
data TechParameters = TechParameters {
    techDelta :: CDouble
    -- ^ Wolfe line search parameter.  Defaults to @0.1@.
    ,techSigma :: CDouble
    -- ^ Wolfe line search parameter.  Defaults to @0.9@.
    ,techGamma :: CDouble
    -- ^ Decay factor for bracket interval width.  Defaults to
    -- @0.66@.
    ,techRho :: CDouble
    -- ^ Growth factor when searching for initial bracketing
    -- interval.  Defaults to @5@.
    ,techEta :: CDouble
    -- ^ Lower bound for the conjugate gradient update parameter
    -- @beta_k@ is @techEta * ||d||_2@.  Defaults to @0.01@.
    ,techPsi0 :: CDouble
    -- ^ Factor used in starting guess for iteration 1.  Defaults
    -- to @0.01@.
    ,techPsi1 :: CDouble
    -- ^ In performing a QuadStep, we evaluate the function at
    -- @psi1 * previous step@.  Defaults to @0.1@.
    ,techPsi2 :: CDouble
    -- ^ When starting a new CG iteration, our initial guess for
    -- the line search stepsize is @psi2 * previous step@.
    -- Defaults to @2@.
    } deriving (Eq, Ord, Show, Read)



-- | How verbose we should be.
data Verbose =
      Quiet
      -- ^ Do not output anything to @stdout@, which most of the
      -- time is good.
    | Verbose
      -- ^ Print what work is being done on each iteraction.
    | VeryVerbose
      -- ^ Print information about every step, may be useful for
      -- troubleshooting.
      deriving (Eq, Ord, Show, Read, Enum)

-- | Line search methods that may be used.
data LineSearch =
      ApproximateWolfe
      -- ^ Use approximate Wolfe line search.
    | AutoSwitch CDouble
      -- ^ Use ordinary Wolfe line search, switch to approximate
      -- Wolfe when
      --
      -- > |f_{k+1} - f_k| < AWolfeFac * C_k
      --
      -- where @C_k@ is the average size of cost and
      -- @AWolfeFac@ is the parameter to this constructor.
      deriving (Eq, Ord, Show, Read)

-- | Stop rules used to decided when to stop iterating.
data StopRules =
      DefaultStopRule CDouble
      -- ^ @DefaultStopRule stop_fac@ stops when
      --
      -- > |g_k|_infty <= max(grad_tol, |g_0|_infty * stop_fac)
      --
      -- where @|g_i|_infty@ is the maximum absolute component of
      -- the gradient at the @i@-th step.
    | AlternativeStopRule
      -- ^ @AlternativeStopRule@ stops when
      --
      -- > |g_k|_infty <= grad_tol * (1 + |f_k|)
      deriving (Eq, Ord, Show, Read)

-- | How to calculate the estimated error in the function value.
data EstimateError =
      AbsoluteEpsilon CDouble
      -- ^ @AbsoluteEpsilon eps@ estimates the error as @eps@.
    | RelativeEpsilon CDouble
      -- ^ @RelativeEpsilon eps@ estimates the error as @eps * C_k@.
      deriving (Eq, Ord, Show, Read)