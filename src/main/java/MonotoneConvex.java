import java.util.Arrays;

public class MonotoneConvex {

    private double[] forwardRates;
    private double[] terms;
    private double[] discreteTerms;
    private double[] discreteForwardRates;  //discrete forward rates
    private double[] values; //rate values from data
    private double[] discreteRatesValues;//discrete forward value
    private double[] discreteInterpolantatNode; //nodes where there are observable rates/fwd rates
    private int indexOfLastTerm; //index of last term
    private boolean Negative_Forwards_Allowed = false;
    private int lastIndexUsed;

    public MonotoneConvex(double[] values, double[] terms) {
        this.values = values;
        this.terms = terms;
        this.indexOfLastTerm = terms.length - 1;
        this.discreteForwardRates = new double[terms.length];
        this.forwardRates = new double[terms.length];
        this.discreteTerms = new double[terms.length]; //discreteTerms(-1 To iIndex)
        this.discreteRatesValues = new double[terms.length];
        this.discreteInterpolantatNode = new double[terms.length];
    }

    public double interpolate(double term) {
        setUpCurves();
        int i;
        if (term <= 0) {
            return this.forwardRates[0];
        } else if (term > this.terms[indexOfLastTerm]) {
            return interpolate(this.terms[this.indexOfLastTerm]) * this.terms[this.indexOfLastTerm] / term + forward(this.terms[this.indexOfLastTerm]) * (1 - this.terms[this.indexOfLastTerm] / term);
        } else {
            i = getIndex(term);
        }
        //the x in (25)
        double x = (term - this.terms[i]) / (this.terms[i + 1] - this.terms[i]);
        double gZero = this.forwardRates[i] - this.discreteForwardRates[i + 1];
        double gOne = this.forwardRates[i + 1] - this.discreteForwardRates[i + 1];
        double gFunction;

        if (x == 0.0 || x == 1.0) {
            gFunction = 0;
        } else if ((gZero < 0 && -0.5 * gZero <= gOne && gOne <= -2 * gZero) || (gZero > 0 && -0.5 * gZero >= gOne && gOne >= -2 * gZero)) {
            //zone (i)
            gFunction = gZero * (x - 2 * Math.pow(x, 2) + Math.pow(x, 3)) + gOne * (-Math.pow(x, 2) + Math.pow(x, 3));
        } else if ((gZero < 0 && gOne > -2 * gZero) || (gZero > 0 && gOne < -2 * gZero)) {
            //zone (ii)
            //(29)
            double eta = (gOne + 2 * gZero) / (gOne - gZero);
            //(28)
            if (x <= eta) {
                gFunction = gZero * x;
            } else {
                gFunction = gZero * x + (gOne - gZero) * Math.pow((x - eta), 3) / Math.pow((1 - eta), 2) / 3;
            }

        } else if ((gZero > 0 && 0 > gOne && gOne > -0.5 * gZero) || (gZero < 0 && 0 < gOne && gOne < -0.5 * gZero)) {
            //zone (iii)
            //(31)
            double eta = 3 * gOne / (gOne - gZero);
            //(30)
            if (x < eta) {
                gFunction = gOne * x - 1.0 / 3.0 * (gZero - gOne) * (Math.pow((eta - x), 3) / Math.pow(eta, 2) - eta);
            } else {
                gFunction = (2.0 / 3.0 * gOne + 1.0 / 3.0 * gZero) * eta + gOne * (x - eta);
            }
        } else if (gZero == 0 || gOne == 0) {
            gFunction = 0;
        } else {
            //zone (iv)
            //(33)
            double eta = gOne / (gOne + gZero);
            //(34)
            double A = -gZero * gOne / (gZero + gOne);
            //(32)
            if (x <= eta) {
                gFunction = A * x - 1.0 / 3.0 * (gZero - A) * (Math.pow((eta - x), 3) / Math.pow(eta, 2) - eta);
            } else {
                gFunction = (2.0 / 3.0 * A + 1.0 / 3.0 * gZero) * eta + A * (x - eta) + (gOne - A) / 3 * Math.pow((x - eta), 3) / Math.pow((1 - eta), 2);
            }
        } //(12)

        return 1 / term * (this.terms[i] * this.discreteInterpolantatNode[i] + (term - this.terms[i]) * this.discreteForwardRates[i + 1] + (this.terms[i + 1] - this.terms[i]) * gFunction);
    }


    public double forward(double Term) {
        this.setUpCurves();
        int i;
        if (Term <= 0) {
            return this.forwardRates[0];
        } else if (Term > this.terms[this.indexOfLastTerm]) {
            return forward(this.terms[this.indexOfLastTerm]);
        } else {
            i = getIndex(Term);
        }
        //the x in (25)
        double x = (Term - terms[i]) / (terms[i + 1] - terms[i]);
        double gZero = forwardRates[i] - discreteForwardRates[i + 1];
        double gOne = forwardRates[i + 1] - discreteForwardRates[i + 1];
        double gFunction;

        if (x == 0) {
            gFunction = gZero;
        } else if (x == 1) {
            gFunction = gOne;
        } else if ((gZero < 0 && -0.5 * gZero <= gOne && gOne <= -2 * gZero) || (gZero > 0 && -0.5 * gZero >= gOne && gOne >= -2 * gZero)) {
            //zone (i)
            gFunction = gZero * (1 - 4 * x + 3 * Math.pow(x, 2)) + gOne * (-2 * x + 3 * Math.pow(x, 2));

        } else if ((gZero < 0 && gOne > -2 * gZero) || (gZero > 0 && gOne < -2 * gZero)) {//zone (ii)
            //(29)
            double eta = (gOne + 2 * gZero) / (gOne - gZero);
            //(28)
            if (x <= eta) {
                gFunction = gZero;
            } else {
                gFunction = gZero + (gOne - gZero) * Math.pow(((x - eta) / (1 - eta)), 2);
            }

        } else if ((gZero > 0 && 0 > gOne && gOne > -0.5 * gZero) || (gZero < 0 && 0 < gOne && gOne < -0.5 * gZero)) {
            //zone (iii)
            //(31)
            double eta = 3 * gOne / (gOne - gZero);
            //(30)
            if (x < eta) {
                gFunction = gOne + (gZero - gOne) * Math.pow(((eta - x) / eta), 2);
            } else {
                gFunction = gOne;
            }
        } else if (gZero == 0 || gOne == 0) {
            gFunction = 0;
        } else {
            //zone (iv)
            //(33)
            double eta = gOne / (gOne + gZero);
            //(34)
            double A = -gZero * gOne / (gZero + gOne);
            //(32)
            if (x <= eta) {
                gFunction = A + (gZero - A) * Math.pow(((eta - x) / eta), 2);
            } else {
                gFunction = A + (gOne - A) * Math.pow(((eta - x) / (1 - eta)), 2);
            }
        }
        //(26)
        return gFunction + this.discreteForwardRates[i + 1];
    }

    public int getIndex(double Term) {
        int iLastIndex = (int) (bound(0, this.lastIndexUsed, this.indexOfLastTerm));

        while (true) {
            if (Term >= this.terms[iLastIndex]) {

                if (iLastIndex == this.terms[this.indexOfLastTerm]) {
                    if (Term == this.terms[iLastIndex]) {
                        this.lastIndexUsed = this.indexOfLastTerm - 1;
                        break;
                    } else {
                        this.lastIndexUsed = this.indexOfLastTerm;
                        break;
                    }
                } else {
                    if (Term >= this.terms[iLastIndex + 1]) {
                        iLastIndex = iLastIndex + 1;
                    } else {
                        this.lastIndexUsed = iLastIndex;
                        break;
                    }
                }
            } else {
                if (iLastIndex == 0) {
                    this.lastIndexUsed = 0;
                    break;
                } else {
                    if (Term >= terms[iLastIndex - 1]) {
                        this.lastIndexUsed = iLastIndex - 1;
                        break;
                    } else {
                        iLastIndex = iLastIndex - 1;
                    }
                }
            }
        }

        return this.lastIndexUsed;
    }

    private void setUpCurves() {
        //extend the curve to time 0, for the purpose of calculating forward at time 1
        this.discreteTerms[0] = 0.0;
        this.discreteRatesValues[0] = discreteRatesValues[1];
        //step 1
        for (int j = 1; j <= this.indexOfLastTerm; j++) {
            this.discreteForwardRates[j] = (this.terms[j] * this.values[j] - this.terms[j - 1] * this.values[j - 1]) / (this.terms[j] - this.terms[j - 1]);
            this.discreteInterpolantatNode[j] = this.values[j];
        }
        for (int j = 1; j <= this.indexOfLastTerm - 1; j++) {
            this.forwardRates[j] = (this.terms[j] - this.terms[j - 1]) / (this.terms[j + 1] - this.terms[j - 1]) * this.discreteForwardRates[j + 1]
                    + (this.terms[j + 1] - this.terms[j]) / (this.terms[j + 1] - this.terms[j - 1]) * this.discreteForwardRates[j];
        }
        //(23)
        this.forwardRates[0] = this.discreteForwardRates[1] - 0.5 * (this.forwardRates[1] - this.discreteForwardRates[1]);
        //(24)
        this.forwardRates[this.indexOfLastTerm] = this.discreteForwardRates[this.indexOfLastTerm] - 0.5 * (this.forwardRates[this.indexOfLastTerm - 1] - this.discreteForwardRates[this.indexOfLastTerm]);
        //step 3
        if (!this.Negative_Forwards_Allowed) {
            this.forwardRates[0] = this.bound(0, this.forwardRates[0], 2 * this.discreteForwardRates[1]);
            for (int j = 1; j <= this.indexOfLastTerm - 1; j++) {
                this.forwardRates[j] = this.bound(0, forwardRates[j], 2 * Math.min(this.discreteForwardRates[j], this.discreteForwardRates[j + 1]));
            }

            this.forwardRates[this.indexOfLastTerm] = this.bound(0, this.forwardRates[this.indexOfLastTerm], 2 * this.discreteForwardRates[this.indexOfLastTerm]);
        } else {
            double termRate = 0;
            for (int j = 1; j < this.indexOfLastTerm - 1; j++) {
                this.discreteForwardRates[j]=this.values[j];
                termRate = termRate + this.discreteForwardRates[j] * (this.terms[j] - this.terms[j - 1]);
                this.discreteInterpolantatNode[j] = termRate / this.terms[j];
            }
        }
    }

    //Helpers
    private double bound(double minimum, double variable, double maximum) {
        if (variable < minimum) {
            return minimum;
        } else if (variable > maximum) {
            return maximum;
        } else {
            return variable;
        }
    }

}
