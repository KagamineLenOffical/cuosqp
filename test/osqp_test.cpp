#include <gflags/gflags.h>
#include <igl/canonical_quaternions.h>
#include <igl/xml/serialize_xml.h>
#include <osqp.h>
#include <types.h>
DEFINE_string(filename_prefix,
              "", "common prefix for an instance of data, default empty");

//void settingsinit(osqp::OsqpSettings &settings) {
//    settings.eps_rel = 1e-4;
//    settings.max_iter = 100000;
//    settings.eps_prim_inf = 1e-6;
//    settings.eps_dual_inf = 1e-6;
//    settings.adaptive_rho = true;
//    settings.warm_start = true;
//}
struct Instance{

};

Eigen::Map<const Eigen::VectorXd> primal_solution(OSQPSolver *solver)  {
    OSQPSolution *solution = solver->solution;
    return Eigen::Map<const Eigen::VectorXd>(solution->x, solver->work->data->n);
}

int main(int argc, char *argv[]) {
    using google::ParseCommandLineFlags;
    ParseCommandLineFlags(&argc, &argv, true);
    Eigen::SparseMatrix<double> constrains, lhs;
    std::string prefix = FLAGS_filename_prefix;
    prefix = R"(D:\linear_solver\cuosqp\cmake-build-release-visual-studio\out\Release\0423_1_test_yrf_same__0_-1_4_)";
    std::cout << "common prefix '"
              << prefix << "'\n";
    igl::xml::deserialize_xml(constrains, "constrains", prefix + "constrain.xml");
    igl::xml::deserialize_xml(lhs, "lhs", prefix + "lhs.xml");
    Eigen::MatrixXd lb, ub, rhs, result;
    igl::xml::deserialize_xml(lb, "lb", prefix + "lb.xml");
    igl::xml::deserialize_xml(ub, "ub", prefix + "ub.xml");
    igl::xml::deserialize_xml(rhs, "rhs", prefix + "rhs.xml");
    igl::xml::deserialize_xml(result, "result", prefix + "result.xml");
    OSQPSettings *settings = (OSQPSettings *)malloc(sizeof(OSQPSettings));
    if(settings)osqp_set_default_settings(settings);
    settings->eps_rel = 1e-4;
    settings->max_iter = 100000;
    settings->eps_prim_inf = 1e-6;
    settings->eps_dual_inf = 1e-6;
    settings->adaptive_rho = true;
    settings->warm_start = true;


    //settings->polish = true;

    Eigen::SparseMatrix<double, Eigen::ColMajor, c_int> objective_matrix;
    Eigen::VectorXd objective_vector;
    Eigen::SparseMatrix<double, Eigen::ColMajor, c_int> constraint_matrix;
    Eigen::VectorXd lower_bounds;
    Eigen::VectorXd upper_bounds;

    objective_matrix = -lhs;
    constraint_matrix = constrains;
    lower_bounds = lb;
    upper_bounds = ub;
    objective_vector = Eigen::VectorXd(rhs);

    Eigen::VectorXd clipped_lower_bounds = lower_bounds.cwiseMax(-OSQP_INFTY);
    Eigen::VectorXd clipped_upper_bounds = upper_bounds.cwiseMin(OSQP_INFTY);
    c_int m = constraint_matrix.rows(), n = constraint_matrix.cols();

    Eigen::SparseMatrix<double, Eigen::ColMajor, c_int>
            objective_matrix_upper_triangle =
            objective_matrix.triangularView<Eigen::Upper>();

    csc *P = reinterpret_cast<csc*>(malloc(sizeof(csc)));
    csc *A = reinterpret_cast<csc*>(malloc(sizeof(csc)));

    csc_set_data(P,n,n,objective_matrix_upper_triangle.outerIndexPtr()[n],
                 const_cast<double*>(objective_matrix_upper_triangle.valuePtr()),
                 const_cast<c_int*>(objective_matrix_upper_triangle.innerIndexPtr()),
                 const_cast<c_int*>(objective_matrix_upper_triangle.outerIndexPtr()));

    csc_set_data(A,m,n,constraint_matrix.outerIndexPtr()[n],
                 const_cast<double*>(constraint_matrix.valuePtr()),
                 const_cast<c_int*>(constraint_matrix.innerIndexPtr()),
                 const_cast<c_int*>(constraint_matrix.outerIndexPtr()));

    OSQPSolver *solver;

    c_float *q = const_cast<double*>(objective_vector.data());
    c_float *l = clipped_lower_bounds.data();
    c_float *u = clipped_upper_bounds.data();
    c_int exitflag = osqp_setup(&solver,P,q,A,l,u,m,n,settings);

    osqp_solve(solver);
    Eigen::MatrixXd optimal_solution = primal_solution(solver);
    std::cout << "Relative difference between pre = "
              << (optimal_solution - result).norm() << std::endl;

    osqp_cleanup(solver);
    free(A);
    free(P);
    free(settings);
    return 0;
    return exitflag;

    //    osqp::OsqpInstance instance;
//    osqp::OsqpSettings settings;
//    osqp::OsqpSolver solver;
//    settingsinit(settings);
//    instance.objective_matrix = -lhs;
//    instance.constraint_matrix = constrains;
//    instance.lower_bounds = lb;
//    instance.upper_bounds = ub;
//    instance.objective_vector = Eigen::VectorXd(rhs);
//    auto status = solver.Init(instance, settings);
//    std::cout << "solver init status "
//              << status << std::endl;
//    osqp::OsqpExitCode exit_code = solver.Solve();
//    std::cout << "solver solve status "
//              << ToString(exit_code) << std::endl;
//    double optimal_objective = solver.objective_value();
//    Eigen::MatrixXd optimal_solution = solver.primal_solution();
//    std::cout << "Relative difference between pre = "
//              << (optimal_solution - result).norm() << std::endl;
}