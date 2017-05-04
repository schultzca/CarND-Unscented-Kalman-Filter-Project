#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
    /**
    * Constructor.
    */
    Tools();

    /**
    * Destructor.
    */
    virtual ~Tools();

    /**
     * String formatter for vectors and arrays that makes it easy to load
     * as numpy array.
     */
    static Eigen::IOFormat NumpyArrayFormat;

    /**
     * A helper method to calculate RMSE.
     */
    Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

};

#endif /* TOOLS_H_ */
