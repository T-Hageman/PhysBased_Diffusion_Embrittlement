function Solve(obj)
	%SOLVE solves a single time increment for the nonlinear system of
	%equations through a Newton-Raphson procedure, combined with a
	%staggered solution scheme

	outerStop = false;
	outerit = 0;
	err0_array = [];

	while outerStop==false		%Staggered scheme convergence loop
		outerit = outerit + 1;
		fprintf("Loop " + string(outerit) + ":\n");

		for stp=1:obj.NSteps	%staggered scheme steps
			fprintf("  SubStep " + string(stp) + "\n");
    		stop = false;
    		it = 0;
		
			stepnum = size(obj.convergence_log, 1)+1;
 		
			%perform once-per-step calculations
			obj.physics.OncePerStep(stp);

			%assemble system matrices
    		obj.physics.Assemble(stp);
    		obj.physics.Constrain(stp);


    		recalc_pre=true;
    		En_err0 = -1;
    		curr_max_it = obj.maxIt;
    		while(stop == false)	%Newton-Raphson solver loop
		
        		fprintf("    Solving it:" + string(it) + "      ");
        		tsolve = tic;
        		
        		recalc_pre = true;
        		if (recalc_pre)
            		%[P,R,C] = equilibrate(obj.physics.K{stp});
            		%recalc_pre = false;
        		end
		
        		if false   %use preconditioned system
            		d = -R*P*obj.physics.fint{stp};
            		B = R*P*obj.physics.K{stp}*C;
					%cond_num(it+1)=condest(B);
					if true
						dy = B\d;
					else
						[L,U] = ilu(B,struct('type','nofill'));
						dy = gmres(B,d,[],1e-4,500,L,U);
					end
            		dx = C*dy;
				else %do not allow any preconditioning
            		dx = -obj.physics.K{stp}\obj.physics.fint{stp};
        		end
        		tsolve = toc(tsolve);
        		fprintf("        (Solver time:"+string(tsolve)+")\n");
				%fprintf("Conditioning numbers: "+string(cond_num(it+1))+"\n");
		
				%line-search
        		if (obj.linesearch && it>-1)
            		e0 = obj.physics.fint{stp}'*dx;
            		obj.physics.Update(dx, stp);
            		
            		obj.physics.Assemble(stp);
            		obj.physics.Constrain(stp);
            		
            		e1 = obj.physics.fint{stp}'*dx;
            		factor = -e0/(e1-e0);
            		factor = max(obj.linesearchLims(1), min(obj.linesearchLims(2), factor));
            		obj.physics.Update(-(1-factor)*dx, stp);
            		fprintf("    Linesearch: " + string(e0) + " -> " + string(e1) + ":  eta=" + string(factor) +"\n");
        		else
            		obj.physics.Update(dx, stp);
        		end
        		
        		% re-assemble system
        		obj.physics.Assemble(stp);
        		obj.physics.Constrain(stp);

				% convergence check
        		if (En_err0 < 0)
            		En_err0 = sum(abs(obj.physics.fint{stp}.*dx));
            		En_err = En_err0;
		
            		if (En_err0==0)
                		En_err0 = 1e-12;
					end

					err0_array(stp) = En_err0;
        		else
            		En_err = sum(abs(obj.physics.fint{stp}.*dx));
        		end
        		En_err_n = En_err/En_err0;
		
				obj.convergence_log(stepnum,stp,it+1) = En_err_n;
                    		
        		fprintf("    Residual:" + string(En_err_n) + "   ("+string(En_err)  +") \n");
        		
        		it=it+1;
        		if (it>curr_max_it || En_err_n<obj.Conv || En_err<obj.tiny)
                	stop = true;
        		end
			end
		end

		if (outerit > obj.OuterLoops || max(abs(err0_array)) < obj.tiny)
			outerStop = true;
		end

	end
    
	%commit time and path-dependent states
	obj.physics.Commit("Pathdep");
    obj.physics.Commit("Timedep");
end

