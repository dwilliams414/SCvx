% Build Mex
codegen -config:mex -v cr3bp_discretization.m -args {coder.typeof(randn(6, 1)), coder.typeof(0.0), coder.typeof(0.2), coder.typeof(randn(3, 1)), coder.typeof(0.012)}