#include "src/audi.hpp"

using namespace audi;

int main() {
	gdual x(0., "y", 3);

	auto f = erf(x);
	std::cout << f << std::endl;
}


.def("__getstate__", [](const gdual &p) {
            // Returns a tuple that contains the string representation of
            // a gdual as obtained from boost serialization
            std::stringstream ss;
            boost::archive::text_oarchive oa(ss);
            oa << p;
            return py::make_tuple(oa.str());
        })
        .def("__setstate__", [](gdual &p, py::tuple t) {
            if (t.size() != 1)
                throw std::runtime_error("Invalid state!");

            // Invoke the default constructor. 
            new (&p) gdual;
            // Reconstruct the gdual
            std::stringstream ss(t[0]);
            boost::archive::text_iarchive ia(ss);
            ia >> p;
        })