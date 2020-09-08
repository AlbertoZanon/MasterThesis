package com.albertozanon.MercurioModel;

import java.time.LocalDate;
import java.time.Month;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;
import java.util.function.DoubleUnaryOperator;

import org.junit.Assert;
import org.junit.Test;

import com.albertozanon.MercurioModel.LIBORCorrelationModelExponentialDecayWithMercurioModification;
import com.albertozanon.MercurioModel.LIBORMarketModelFromCovarianceModelWithMercurioModification;
import com.albertozanon.MercurioModel.LIBORVolatilityModelFourParameterExponentialFormWithMercurioModification;
import com.albertozanon.MercurioModel.VolatilityReductionMercurioModel;

import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.calibration.ParameterObject;
import net.finmath.marketdata.calibration.Solver;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelFromCurvesAndVols;
import net.finmath.marketdata.model.curves.Curve;
import net.finmath.marketdata.model.curves.DiscountCurve;
import net.finmath.marketdata.model.curves.DiscountCurveFromForwardCurve;
import net.finmath.marketdata.model.curves.DiscountCurveInterpolation;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.marketdata.model.curves.ForwardCurveFromDiscountCurve;
import net.finmath.marketdata.model.curves.CurveInterpolation.ExtrapolationMethod;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationEntity;
import net.finmath.marketdata.model.curves.CurveInterpolation.InterpolationMethod;
import net.finmath.marketdata.products.AnalyticProduct;
import net.finmath.marketdata.products.Swap;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionFromMersenneRandomNumbers;
import net.finmath.montecarlo.RandomVariableFactory;
import net.finmath.montecarlo.RandomVariableFromArrayFactory;
import net.finmath.montecarlo.assetderivativevaluation.MonteCarloAssetModel;
import net.finmath.montecarlo.assetderivativevaluation.models.BlackScholesModel;
import net.finmath.montecarlo.assetderivativevaluation.products.EuropeanOption;
import net.finmath.montecarlo.interestrate.LIBORMarketModel;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationModel;
import net.finmath.montecarlo.interestrate.LIBORMonteCarloSimulationFromLIBORModel;
import net.finmath.montecarlo.interestrate.models.covariance.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.models.covariance.DisplacedLocalVolatilityModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCorrelationModel;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORCovarianceModelFromVolatilityAndCorrelation;
import net.finmath.montecarlo.interestrate.models.covariance.LIBORVolatilityModel;
import net.finmath.montecarlo.process.EulerSchemeFromProcessModel;
import net.finmath.optimizer.SolverException;
import net.finmath.plots.Named;
import net.finmath.plots.Plot;
import net.finmath.plots.Plot2D;
import net.finmath.plots.PlotableFunction2D;
import net.finmath.plots.Plots;
import net.finmath.stochastic.RandomVariable;
import net.finmath.time.Schedule;
import net.finmath.time.ScheduleGenerator;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationFromArray;
import net.finmath.time.businessdaycalendar.BusinessdayCalendarExcludingTARGETHolidays;

/**
 * Testing Plot2D.
 *
 * @author Chrisitan Fries
 */
public class PlotPaths {
	@Test
	public void test2() throws CalculationException, InterruptedException {


		final double modelInitialValue = 100.0;
		final double modelRiskFreeRate = 0.05;
		final double modelVolatility = 0.20;

		final double maturity = 3.0;
		final double strike = 106.0;

		// Create a model
		final var model = new BlackScholesModel(modelInitialValue, modelRiskFreeRate, modelVolatility);

		// Create a corresponding MC process
		final var td = new TimeDiscretizationFromArray(0.0, 300, 0.01);
		final var brownianMotion = new BrownianMotionFromMersenneRandomNumbers(td, 1, 100000, 3231);
		final var process = new EulerSchemeFromProcessModel(model, brownianMotion);

		// Using the process (Euler scheme), create an MC simulation of a Black-Scholes model
		final var simulation = new MonteCarloAssetModel(process);

		final EuropeanOption europeanOption = new EuropeanOption(maturity, strike);

		try {
			final Plot plot = Plots.createHistogramBehindValues(simulation.getAssetValue(maturity, 0 /* assetIndex */), europeanOption.getValue(0.0, simulation), 100, 5.0);
			plot.show();
			Thread.sleep(20000);
		} catch (final Exception e) {
			e.printStackTrace();
		}
	}


	@Test
	public void test() throws SolverException, CalculationException {
		RandomVariableFactory randomVariableFactory = new RandomVariableFromArrayFactory();
		final int numberOfPaths		= 10; // ideale è 10'000 paths!
		final int numberOfFactors	= 1;

		final long millisCurvesStart = System.currentTimeMillis();


		System.out.println("Calibration to Swaptions.\n");

		System.out.println("Calibration of rate curves:");

		final AnalyticModel curveModel = getCalibratedCurve();

		// Create the forward curve (initial value of the LIBOR market model)
		final ForwardCurve forwardCurve = curveModel.getForwardCurve("ForwardCurveFromDiscountCurve(discountCurve-EUR,1D)");
		final DiscountCurve discountCurve = curveModel.getDiscountCurve("discountCurve-EUR");
		//		curveModel.addCurve(discountCurve.getName(), discountCurve);
		final double lastTime	= 1.0;
		final double dtLibor	= 0.2;
		final double dt	= 0.001;
		final TimeDiscretization timeDiscretizationFromArray = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dt), dt);
		final TimeDiscretization liborPeriodDiscretization = new TimeDiscretizationFromArray(0.0, (int) (lastTime / dtLibor), dtLibor);
		final BrownianMotion brownianMotion = new net.finmath.montecarlo.BrownianMotionLazyInit(timeDiscretizationFromArray, numberOfFactors, numberOfPaths, 2109 /* seed */); 

		// Create a model
		//	final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTimeHomogenousPiecewiseConstantWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, timeToMaturityDiscretization, arrayValues);
		LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelFourParameterExponentialFormWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, 0.002, 0.005, 0.2, 0.002, true); //0.20/100.0, 0.05/100.0, 0.10, 0.05/100.0, 
		//final LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelPiecewiseConstantWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization,optionMaturityDiscretization,timeToMaturityDiscretization, 0.50 / 100);
		
		final LIBORCorrelationModel correlationModel = new LIBORCorrelationModelExponentialDecayWithMercurioModification(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, 0.05, false);
		//AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelExponentialForm5Param(timeDiscretizationFromArray, liborPeriodDiscretization, numberOfFactors, new double[] { 0.20/100.0, 0.05/100.0, 0.10, 0.05/100.0, 0.10} );
		final AbstractLIBORCovarianceModelParametric covarianceModelParametric = new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretizationFromArray, liborPeriodDiscretization, volatilityModel, correlationModel);

		// Create blended local volatility model with fixed parameter (0=lognormal, > 1 = almost a normal model).
		final AbstractLIBORCovarianceModelParametric covarianceModelDisplaced = new DisplacedLocalVolatilityModel(covarianceModelParametric, 1.0/0.25, false /* isCalibrateable */);
		final AbstractLIBORCovarianceModelParametric covarianceModelReducedVolatility = new VolatilityReductionMercurioModel(covarianceModelDisplaced);
		
		Map<String, String> properties2 = new HashMap<String, String >();
		properties2.put("measure", LIBORMarketModelFromCovarianceModelWithMercurioModification.Measure.SPOT.name());
		properties2.put("stateSpace", LIBORMarketModelFromCovarianceModelWithMercurioModification.StateSpace.NORMAL.name());

		final LIBORMarketModel MercurioModelNONcalibrated=  new LIBORMarketModelFromCovarianceModelWithMercurioModification(
				liborPeriodDiscretization,
				curveModel,
				forwardCurve,
				new DiscountCurveFromForwardCurve(forwardCurve),
				randomVariableFactory,
				covarianceModelReducedVolatility,
				properties2
				);

		final EulerSchemeFromProcessModel process2 = new EulerSchemeFromProcessModel(MercurioModelNONcalibrated,brownianMotion);
		final LIBORModelMonteCarloSimulationModel simulationMercurioModelNONcalibrated = new LIBORMonteCarloSimulationFromLIBORModel(process2);
	

		final DoubleUnaryOperator function = (time) -> {


			RandomVariable z = null;
			try {
				z = simulationMercurioModelNONcalibrated.getLIBOR(time,0.75,1.0);
			} catch (CalculationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			double path = z.get(7);
			return path;
		};

		final DoubleUnaryOperator function2 = (time) -> {


			RandomVariable z = null;
			try {
				z = simulationMercurioModelNONcalibrated.getLIBOR(time,0.75,1.0);
			} catch (CalculationException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			double path = z.get(5);
			return path;
		};
		
		Plot plot = new Plot2D(List.of(new PlotableFunction2D(0.0, 1.0, 1000, function), new PlotableFunction2D(0.0, 1.0, 1000, function2)));
		//final Plot plot = new Plot2D(0.0, 1.0, 1000, Arrays.asList(
				//new Named<DoubleUnaryOperator>("Maturity 1", function)));

		plot.setTitle("MC Simulation").setXAxisLabel("Time").setYAxisLabel("Rate's Value").setIsLegendVisible(true);
		try {
			plot.show();
			Thread.sleep(20000);

		} catch (final Exception e) {
		}
	}
	public AnalyticModel getCalibratedCurve() throws SolverException {
		final String[] maturity					= { "7D", "14D", "21D", "1M", "2M", "3M", "4M", "5M", "6M", "7M", "8M", "9M", "12M", "15M", "18M", "21M", "2Y", "3Y", "4Y", "5Y", "6Y", "7Y", "8Y", "9Y", "10Y", "11Y", "12Y", "15Y", "20Y", "25Y", "30Y", "40Y", "50Y" };
		final String[] frequency				= { "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
		final String[] frequencyFloat			= { "tenor", "tenor",  "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "tenor", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual", "annual"};
		final String[] daycountConventions	    = { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
		final String[] daycountConventionsFloat	= { "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360", "ACT/360"};
		final double[] rates					= { -0.00553, -0.00553, -0.00553, -0.00553, -0.00555, -0.00556, -0.00559, -0.00564, -0.00568, -0.00572, -0.00577, -0.00581, -0.00592, -0.00601, -0.00608, -0.00613, -0.00619, -0.00627, -0.00622, -0.00606, -0.00582, -0.00553, -0.00519, -0.00482, -0.00442, -0.00402, -0.00362, -0.00261, -0.00189, -0.00197, -0.0023, -0.00286, -0.00333};
		final HashMap<String, Object> parameters = new HashMap<>();

		parameters.put("referenceDate", LocalDate.of(2020, Month.JULY, 31));
		parameters.put("currency", "EUR");
		parameters.put("forwardCurveTenor", "1D");
		parameters.put("maturities", maturity);
		parameters.put("fixLegFrequencies", frequency);
		parameters.put("floatLegFrequencies", frequencyFloat);
		parameters.put("fixLegDaycountConventions", daycountConventions);
		parameters.put("floatLegDaycountConventions", daycountConventionsFloat);
		parameters.put("rates", rates);

		return getCalibratedCurve(null, parameters);
	}

	private static AnalyticModel getCalibratedCurve(final AnalyticModel model2, final Map<String, Object> parameters) throws SolverException {

		final LocalDate	referenceDate		= (LocalDate) parameters.get("referenceDate");
		final String	currency			= (String) parameters.get("currency");
		final String	forwardCurveTenor	= (String) parameters.get("forwardCurveTenor");
		final String[]	maturities			= (String[]) parameters.get("maturities");
		final String[]	frequency			= (String[]) parameters.get("fixLegFrequencies");
		final String[]	frequencyFloat		= (String[]) parameters.get("floatLegFrequencies");
		final String[]	daycountConventions	= (String[]) parameters.get("fixLegDaycountConventions");
		final String[]	daycountConventionsFloat	= (String[]) parameters.get("floatLegDaycountConventions");
		final double[]	rates						= (double[]) parameters.get("rates");

		Assert.assertEquals(maturities.length, frequency.length);
		Assert.assertEquals(maturities.length, daycountConventions.length);
		Assert.assertEquals(maturities.length, rates.length);

		Assert.assertEquals(frequency.length, frequencyFloat.length);
		Assert.assertEquals(daycountConventions.length, daycountConventionsFloat.length);

		final int		spotOffsetDays = 2;
		final String	forwardStartPeriod = "0D";

		final String curveNameDiscount = "discountCurve-" + currency;

		/*
		 * We create a forward curve by referencing the same discount curve, since
		 * this is a single curve setup.
		 *
		 * Note that using an independent NSS forward curve with its own NSS parameters
		 * would result in a problem where both, the forward curve and the discount curve
		 * have free parameters.
		 */
		final ForwardCurve forwardCurve		= new ForwardCurveFromDiscountCurve(curveNameDiscount, referenceDate, forwardCurveTenor);

		// Create a collection of objective functions (calibration products)
		final Vector<AnalyticProduct> calibrationProducts = new Vector<>();
		final double[] curveMaturities	= new double[rates.length+1];
		final double[] curveValue			= new double[rates.length+1];
		final boolean[] curveIsParameter	= new boolean[rates.length+1];
		curveMaturities[0] = 0.0;
		curveValue[0] = 1.0;
		curveIsParameter[0] = false;
		for(int i=0; i<rates.length; i++) {

			final Schedule schedulePay = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequency[i], daycountConventions[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);
			final Schedule scheduleRec = ScheduleGenerator.createScheduleFromConventions(referenceDate, spotOffsetDays, forwardStartPeriod, maturities[i], frequencyFloat[i], daycountConventionsFloat[i], "first", "following", new BusinessdayCalendarExcludingTARGETHolidays(), -2, 0);

			curveMaturities[i+1] = Math.max(schedulePay.getPayment(schedulePay.getNumberOfPeriods()-1),scheduleRec.getPayment(scheduleRec.getNumberOfPeriods()-1));
			curveValue[i+1] = 1.0;
			curveIsParameter[i+1] = true;
			calibrationProducts.add(new Swap(schedulePay, null, rates[i], curveNameDiscount, scheduleRec, forwardCurve.getName(), 0.0, curveNameDiscount));
		}

		final InterpolationMethod interpolationMethod = InterpolationMethod.LINEAR;

		// Create a discount curve
		final DiscountCurveInterpolation discountCurveInterpolation = DiscountCurveInterpolation.createDiscountCurveFromDiscountFactors(
				curveNameDiscount								/* name */,
				referenceDate	/* referenceDate */,
				curveMaturities	/* maturities */,
				curveValue		/* discount factors */,
				curveIsParameter,
				interpolationMethod ,
				ExtrapolationMethod.CONSTANT,
				InterpolationEntity.LOG_OF_VALUE
				);

		/*
		 * Model consists of the two curves, but only one of them provides free parameters.
		 */
		AnalyticModel model = new AnalyticModelFromCurvesAndVols(new Curve[] { discountCurveInterpolation, forwardCurve });

		/*
		 * Create a collection of curves to calibrate
		 */
		final Set<ParameterObject> curvesToCalibrate = new HashSet<>();
		curvesToCalibrate.add(discountCurveInterpolation);

		/*
		 * Calibrate the curve
		 */
		final Solver solver = new Solver(model, calibrationProducts, 0.0, 1E-4 /* target accuracy */);
		final AnalyticModel calibratedModel = solver.getCalibratedModel(curvesToCalibrate);
		System.out.println("Solver reported acccurary....: " + solver.getAccuracy());

		Assert.assertEquals("Calibration accurarcy", 0.0, solver.getAccuracy(), 1E-3);

		// Get best parameters
		final double[] parametersBest = calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getParameter();

		// Test calibration
		model			= calibratedModel;

		double squaredErrorSum = 0.0;
		for(final AnalyticProduct c : calibrationProducts) {
			final double value = c.getValue(0.0, model);
			final double valueTaget = 0.0;
			final double error = value - valueTaget;
			squaredErrorSum += error*error;
		}
		final double rms = Math.sqrt(squaredErrorSum/calibrationProducts.size());

		System.out.println("Independent checked acccurary: " + rms);

		System.out.println("Calibrated discount curve: ");
		for(int i=0; i<curveMaturities.length; i++) {
			final double maturity = curveMaturities[i];
			System.out.println(maturity + "\t" + calibratedModel.getDiscountCurve(discountCurveInterpolation.getName()).getDiscountFactor(maturity));
		}
		return model;
	}
}