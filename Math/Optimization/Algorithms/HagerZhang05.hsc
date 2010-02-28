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

    ,qdecay :: Double
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

    ,quadraticStep :: Maybe Double
    -- ^ When to attempt quadratic interpolation in line search.
    -- If @Nothing@ then never try a quadratic interpolation
    -- step.  If @Just cutoff@, then attemp quadratic
    -- interpolation in line search when @|f_{k+1} - f_k| / f_k
    -- <= cutoff@.  Defaults to @Just 1e-12@.

    ,debugTol :: Maybe Double
    -- ^ If @Just tol@, then always check that @f_{k+1} - f_k <=
    -- tol * C_k@. Otherwise, if @Nothing@ then no checking of
    -- function values is done.  Defaults to @Nothing@.

    ,initialStep :: Maybe Double
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

    ,funcEpsilon :: Double
    -- ^ Stop when @-alpha * dphi0@, the estimated change in
    -- function value, is less than @funcEpsilon * |f|@.
    -- Defaults to @0@.

    ,nanRho :: Double
    -- ^ After encountering @NaN@ while calculating the step
    -- length, growth factor when searching for a bracketing
    -- interval.  Defaults to @1.3@.

    ,techParameters :: TechParameters
    -- ^ Technical parameters which you probably should not
    -- touch.
    } deriving (Eq, Ord, Show, Read)

instance Storable Parameters where
    sizeOf _    = #{size cg_parameter}
    alignment _ = alignment (undefined :: Double)
    peek ptr    = do
      v_printFinal    <- #{peek cg_parameter, PrintFinal}  ptr
      v_printParams   <- #{peek cg_parameter, PrintParms}  ptr
      v_verbose       <- #{peek cg_parameter, PrintLevel}  ptr
      v_awolfe        <- #{peek cg_parameter, AWolfe}      ptr
      v_awolfefac     <- #{peek cg_parameter, AWolfeFac}   ptr
      v_qdecay        <- #{peek cg_parameter, Qdecay}      ptr
      v_stopRule      <- #{peek cg_parameter, StopRule}    ptr
      v_stopRuleFac   <- #{peek cg_parameter, StopFac}     ptr
      v_estimateError <- #{peek cg_parameter, PertRule}    ptr
      v_estimateEps   <- #{peek cg_parameter, eps}         ptr
      v_quadraticStep <- #{peek cg_parameter, QuadStep}    ptr
      v_quadraticCut  <- #{peek cg_parameter, QuadCutOff}  ptr
      v_debug         <- #{peek cg_parameter, debug}       ptr
      v_debugTol      <- #{peek cg_parameter, debugtol}    ptr
      v_initialStep   <- #{peek cg_parameter, step}        ptr
      v_maxItersFac   <- #{peek cg_parameter, maxit_fac}   ptr
      v_nexpand       <- #{peek cg_parameter, nexpand}     ptr
      v_nsecant       <- #{peek cg_parameter, nsecant}     ptr
      v_restartFac    <- #{peek cg_parameter, restart_fac} ptr
      v_funcEpsilon   <- #{peek cg_parameter, feps}        ptr
      v_nanRho        <- #{peek cg_parameter, nan_rho}     ptr

      v_delta         <- #{peek cg_parameter, delta}       ptr
      v_sigma         <- #{peek cg_parameter, sigma}       ptr
      v_gamma         <- #{peek cg_parameter, gamma}       ptr
      v_rho           <- #{peek cg_parameter, rho}         ptr
      v_eta           <- #{peek cg_parameter, eta}         ptr
      v_psi0          <- #{peek cg_parameter, psi0}        ptr
      v_psi1          <- #{peek cg_parameter, psi1}        ptr
      v_psi2          <- #{peek cg_parameter, psi2}        ptr

      let tech = TechParameters {techDelta = v_delta
                                ,techSigma = v_sigma
                                ,techGamma = v_gamma
                                ,techRho   = v_rho
                                ,techEta   = v_eta
                                ,techPsi0  = v_psi0
                                ,techPsi1  = v_psi1
                                ,techPsi2  = v_psi2}

      let b :: CInt -> Bool; b = (/= 0)

      return Parameters {printFinal     = b v_printFinal
                        ,printParams    = b v_printParams
                        ,verbose        = case v_verbose :: CInt of
                                            0 -> Quiet
                                            1 -> Verbose
                                            _ -> VeryVerbose
                        ,lineSearch     = if b v_awolfe
                                          then ApproximateWolfe
                                          else AutoSwitch v_awolfefac
                        ,qdecay         = v_qdecay
                        ,stopRules      = if b v_stopRule
                                          then DefaultStopRule v_stopRuleFac
                                          else AlternativeStopRule
                        ,estimateError  = if b v_estimateError
                                          then RelativeEpsilon v_estimateEps
                                          else AbsoluteEpsilon v_estimateEps
                        ,quadraticStep  = if b v_quadraticStep
                                          then Just v_quadraticCut
                                          else Nothing
                        ,debugTol       = if b v_debug
                                          then Just v_debugTol
                                          else Nothing
                        ,initialStep    = case v_initialStep of
                                            0 -> Nothing
                                            x -> Just x
                        ,maxItersFac    = v_maxItersFac
                        ,nexpand        = v_nexpand
                        ,nsecant        = v_nsecant
                        ,restartFac     = v_restartFac
                        ,funcEpsilon    = v_funcEpsilon
                        ,nanRho         = v_nanRho
                        ,techParameters = tech}
    poke ptr p = do
      let i b = if b p then 1 else (0 :: CInt)
          m b = maybe (0 :: CInt) (const 1) (b p)
      #{poke cg_parameter, PrintFinal}  ptr (i printFinal)
      #{poke cg_parameter, PrintParms}  ptr (i printParams)
      #{poke cg_parameter, PrintLevel}  ptr (case verbose p of
                                               Quiet       -> 0 :: CInt
                                               Verbose     -> 1
                                               VeryVerbose -> 3)
      let (awolfe, awolfefac) = case lineSearch p of
                                  ApproximateWolfe -> (1, 0)
                                  AutoSwitch x     -> (0, x)
      #{poke cg_parameter, AWolfe}      ptr (awolfe :: CInt)
      #{poke cg_parameter, AWolfeFac}   ptr awolfefac
      #{poke cg_parameter, Qdecay}      ptr (qdecay p)
      let (stopRule, stopRuleFac) = case stopRules p of
                                      DefaultStopRule x   -> (1, x)
                                      AlternativeStopRule -> (0, 0)
      #{poke cg_parameter, StopRule}    ptr (stopRule :: CInt)
      #{poke cg_parameter, StopFac}     ptr stopRuleFac
      let (pertRule, eps) = case estimateError p of
                              RelativeEpsilon x -> (1,x)
                              AbsoluteEpsilon x -> (0,x)
      #{poke cg_parameter, PertRule}    ptr (pertRule :: CInt)
      #{poke cg_parameter, eps}         ptr eps
      #{poke cg_parameter, QuadStep}    ptr (m quadraticStep)
      #{poke cg_parameter, QuadCutOff}  ptr (maybe 0 id $ quadraticStep p)
      #{poke cg_parameter, debug}       ptr (m debugTol)
      #{poke cg_parameter, debugtol}    ptr (maybe 0 id $ debugTol p)
      #{poke cg_parameter, step}        ptr (maybe 0 id $ initialStep p)
      #{poke cg_parameter, maxit_fac}   ptr (maxItersFac p)
      #{poke cg_parameter, nexpand}     ptr (nexpand p)
      #{poke cg_parameter, nsecant}     ptr (nsecant p)
      #{poke cg_parameter, restart_fac} ptr (restartFac p)
      #{poke cg_parameter, feps}        ptr (funcEpsilon p)
      #{poke cg_parameter, nan_rho}     ptr (nanRho p)

      #{poke cg_parameter, delta}       ptr (techDelta $ techParameters p)
      #{poke cg_parameter, sigma}       ptr (techSigma $ techParameters p)
      #{poke cg_parameter, gamma}       ptr (techGamma $ techParameters p)
      #{poke cg_parameter, rho}         ptr (techRho   $ techParameters p)
      #{poke cg_parameter, eta}         ptr (techEta   $ techParameters p)
      #{poke cg_parameter, psi0}        ptr (techPsi0  $ techParameters p)
      #{poke cg_parameter, psi1}        ptr (techPsi1  $ techParameters p)
      #{poke cg_parameter, psi2}        ptr (techPsi2  $ techParameters p)




-- | Technical parameters which you probably should not touch.
-- You should read the papers of @CG_DESCENT@ to understand how
-- you can tune these parameters.
data TechParameters = TechParameters {
    techDelta :: Double
    -- ^ Wolfe line search parameter.  Defaults to @0.1@.
    ,techSigma :: Double
    -- ^ Wolfe line search parameter.  Defaults to @0.9@.
    ,techGamma :: Double
    -- ^ Decay factor for bracket interval width.  Defaults to
    -- @0.66@.
    ,techRho :: Double
    -- ^ Growth factor when searching for initial bracketing
    -- interval.  Defaults to @5@.
    ,techEta :: Double
    -- ^ Lower bound for the conjugate gradient update parameter
    -- @beta_k@ is @techEta * ||d||_2@.  Defaults to @0.01@.
    ,techPsi0 :: Double
    -- ^ Factor used in starting guess for iteration 1.  Defaults
    -- to @0.01@.
    ,techPsi1 :: Double
    -- ^ In performing a QuadStep, we evaluate the function at
    -- @psi1 * previous step@.  Defaults to @0.1@.
    ,techPsi2 :: Double
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
      DefaultStopRule Double
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
      AbsoluteEpsilon Double
      -- ^ @AbsoluteEpsilon eps@ estimates the error as @eps@.
    | RelativeEpsilon Double
      -- ^ @RelativeEpsilon eps@ estimates the error as @eps * C_k@.
      deriving (Eq, Ord, Show, Read)