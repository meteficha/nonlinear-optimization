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
    -- ^ Print final statistics to @stdout@.

    ,printParams :: Bool
    -- ^ Print parameters to @stdout@ before starting.

    ,verbose :: Verbose
    -- ^ How verbose we should be while computing.  Everything is
    -- printed to @stdout@.

    ,lineSearch :: LineSearch
    -- ^ What kind of line search should be used.

    ,qdecay :: Double
    -- ^ Factor in @[0, 1]@ used to compute average cost
    -- magnitude @C_k@ as follows:
    --
    -- > Q_k = 1 + (qdecay)Q_{k-1},   Q_0 = 0
    -- > C_k = C_{k-1} + (|f_k| - C_{k-1})/Q_k

    ,stopRules :: StopRules
    -- ^ Stop rules that define when the iterations should end.

    ,estimateError :: EstimateError
    -- ^ How to calculate the estimated error in the function value.

    ,quadraticStep :: Maybe Double
    -- ^ When to attempt quadratic interpolation in line search.
    -- If @Nothing@ then never try a quadratic interpolation
    -- step.  If @Just cutoff@, then attemp quadratic
    -- interpolation in line search when @|f_{k+1} - f_k| / f_k
    -- <= cutoff@.

    ,debugTol :: Maybe Double
    -- ^ If @Just tol@, then always check that @f_{k+1} - f_k <=
    -- tol * C_k@. Otherwise, if @Nothing@ then no checking of
    -- function values is done.

    ,initialStep :: Maybe Double
    -- ^ If @Just step@, then use @step@ as the initial step of
    -- the line search.  Otherwise, if @Nothing@ then the initial
    -- step is programatically calculated.

    ,maxItersFac :: Int
    -- ^ Defines the maximum number of iterations.  The process
    -- is aborted when @maxItersFac * n@ iterations are done, where
    -- @n@ is the number of dimensions.

    ,nexpand :: Int
    -- ^ Maximum number of times the bracketing interval grows or
    -- shrinks in the line search.

    ,nsecant :: Int
    -- ^ Maximum number of secant iterations in line search.

    ,restartFac :: Int
    -- ^ Restart the conjugate gradient method after @restartFac
    -- * n@ iterations.

    ,funcEpsilon :: Double
    -- ^ Stop when @-alpha * dphi0@, the estimated change in
    -- function value, is less than @funcEpsilon * |f|@.

    ,nanRho :: Double
    -- ^ After encountering @NaN@ while calculating the step
    -- length, growth factor when searching for a bracketing
    -- interval.

    ,techParameters :: TechParameters
    -- ^ Technical parameters which you probably should not
    -- touch.
    } deriving (Eq, Ord, Show, Read)

-- | Technical parameters which you probably should not touch.
-- You should read the papers of @CG_DESCENT@ to understand how
-- you can tune these parameters.
data TechParameters = TechParameters {
    techDelta :: Double
    -- ^ Wolfe line search parameter.
    ,techSigma :: Double
    -- ^ Wolfe line search parameter.
    ,techGamma :: Double
    -- ^ Decay factor for bracket interval width.
    ,techRho :: Double
    -- ^ Growth factor when searching for initial bracketing
    -- interval.
    ,techEta :: Double
    -- ^ Lower bound for the conjugate gradient update parameter
    -- @beta_k@ is @techEta * ||d||_2@.
    ,techPsi0 :: Double
    -- ^ Factor used in starting guess for iteration 1.
    ,techPsi1 :: Double
    -- ^ In performing a QuadStep, we evaluate the function at
    -- @psi1 * previous step@.
    ,techPsi2 :: Double
    -- ^ When starting a new CG iteration, our initial guess for
    -- the line search stepsize is @psi2 * previous step@.
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
    | AutoSwitch Double
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
      DefaultStopRule Double Double
      -- ^ @DefaultStopRule grad_tol stop_fac@ stops when
      --
      -- > |g_k|_infty <= max(grad_tol, |g_0|_infty * stop_fac)
      --
      -- where @|g_i|_infty@ is the maximum absolute component of
      -- the gradient at the @i@-th step.
    | AlternativeStopRule Double
      -- ^ @AlternativeStopRule grad_tol@ stops when
      --
      -- > |g_k|_infty <= grad_tol * (1 + |f_k|)
      deriving (Eq, Ord, Show, Read)

-- | How to calculate the estimated error in the function value.
data EstimateError =
      AbsoluteEpsilon Double
      -- ^ @AbsoluteEpsilon eps@ estimates the error as @eps@.
    | RelativeEpsilon Double
      -- ^ @RelativeEpsilon eps@ estimates the error as @eps * C_k@.
      deriving (Eq, Ord, Show, Read)