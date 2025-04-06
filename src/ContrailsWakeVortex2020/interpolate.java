package ContrailsWakeVortex2020;

public class interpolate {

    public static double linear(double[] x, double[] y, double x_neu) {
        return linear(x, y, x_neu, false);
    }

    public static double cubicspline(double[] x, double[] y, double x_neu, boolean extrapolate) {
        double  interpolatedvalue=Double.NaN;
        try{
            org.apache.commons.math3.analysis.interpolation.SplineInterpolator splineInterpolator=new org.apache.commons.math3.analysis.interpolation.SplineInterpolator();
            interpolatedvalue = splineInterpolator.interpolate(x, y).value(x_neu);
        }catch (Exception e){
            //System.out.println("cubic interpolation/extrapolation impossible->jump to linear interpolation (value:" + x_neu+", lower bound: "+x[0] +", upper bound: "+x[x.length-1]+")");
            interpolatedvalue=linear(x,y,x_neu,extrapolate);
        }

        return interpolatedvalue;
    }

    public static double linear_windheading(double[] x, double[] y, double x_neu, boolean extrapolate) {
        int n = 0;
        if (x[x.length - 1] <= x_neu) {
            if (!extrapolate) {
                return y[y.length - 1];
            }
            n = x.length - 1;
            // set params for extrapolate above
        } else if (x[0] >= x_neu) {
            if (!extrapolate) {
                return y[0];
            }
            n = 1;
            // set params for extrapolate below
        } else {
            while (x[n] < x_neu) {
                ++n;
            }
        }
        double h2 = x[n];
        double h1 = x[n - 1];
        double f2 = y[n];
        double f1 = y[n - 1];
        //h2 always greater than h1

        double amount = 1 - (x_neu - h1) / (h2 - h1);
        double result = lerpDegrees(f2, f1, amount);
        return result;
        //return f1+((f2-f1)/(h2-h1))*(x_neu-h1);
    }

    public static double lerpDegrees(double start, double end, double amount) {
        double difference = Math.abs(end - start);
        if (difference > 180) {
            // We need to add on to one of the values.
            if (end > start) {
                // We'll add it on to start...
                start += 360;
            } else {
                // Add it on to end.
                end += 360;
            }
        }

        // Interpolate it.
        double value = (start + ((end - start) * amount));

        // Wrap it..
        float rangeZero = 360;

        if (value >= 0 && value <= 360)
            return value;

        return (value % rangeZero);
    }

    public static double linear(double[] x, double[] y, double x_neu, boolean extrapolate) {
        int n = 0;
        if (x[x.length - 1] <= x_neu) {
            if (!extrapolate) {
                return y[y.length - 1];
            }
            n = x.length - 1;
            // set params for extrapolate above
        } else if (x[0] >= x_neu) {
            if (!extrapolate) {
                return y[0];
            }
            n = 1;
            // set params for extrapolate below
        } else {
            while (x[n] < x_neu) {
                ++n;
            }
        }

        double h2 = x[n];
        double h1 = x[n - 1];
        double f2 = y[n];
        double f1 = y[n - 1];

        return f1+((f2-f1)/(h2-h1))*(x_neu-h1);
        //double a = (f2 - f1) / (h2 - h1);
        //return f1 + (x_neu - h1) * a;
    }

    public static double loglog(double[] x, double[] y, double x_neu) {
        // interpoliert linear im loglog Diagramm
        int n = 0;
        if (x_neu > x[x.length - 1]) {
            return y[y.length - 1]; // keine extrapolation nach oben
        }
        while (x[n] < x_neu && n < x.length - 1) {
            n = n + 1;
        }
        if (n == 0) {
            return y[0];
        }
        double h2 = Math.log(x[n]);
        double h1 = Math.log(x[n - 1]);
        double f2 = Math.log(y[n]);
        double f1 = Math.log(y[n - 1]);

        double a = (f2 - f1) / (h2 - h1);
        return Math.exp(f1 + (Math.log(x_neu) - h1) * a);
    }
}
