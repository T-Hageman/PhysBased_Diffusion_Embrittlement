function Solve(obj)
	outerStop = false;
	outerit = 0;
	err0_array = [];

	while outerStop==false
		outerit = outerit + 1;
		fprintf("Loop " + string(outerit) + ":\n");

		for stp=1:obj.NSteps
			fprintf("  SubStep " + string(stp) + "\n");
    		stop = false;
    		it = 0;
		
			stepnum = size(obj.convergence_log, 1)+1;
			%cf = figure(9731);
 		
			obj.physics.OncePerStep(stp);
    		obj.physics.Assemble(stp);
    		obj.physics.Constrain(stp);
    		recalc_pre=true;
    		En_err0 = -1;
    		curr_max_it = obj.maxIt;
    		while(stop == false)
		
        		fprintf("    Solving it:" + string(it) + "      ");
        		tsolve = tic;
        		
        		recalc_pre = true;
        		if (recalc_pre)
            		%[P,R,C] = equilibrate(obj.physics.K{stp});
            		%recalc_pre = false;
        		end
		
        		if false
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
        		else
            		dx = -obj.physics.K{stp}\obj.physics.fint{stp};
        		end
        		tsolve = toc(tsolve);
        		fprintf("        (Solver time:"+string(tsolve)+")\n");
				%fprintf("Conditioning numbers: "+string(cond_num(it+1))+"\n");
		
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
        		
        		% convergence
        		obj.physics.Assemble(stp);
        		obj.physics.Constrain(stp);
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
            		%obj.physics.Commit("Pathdep");
            		%irr = obj.physics.Irreversibles();
            		%if (irr == false)
                		stop = true;
	%             	else
	%                 	obj.physics.Assemble();
	%                 	obj.physics.Constrain();
	%                 	recalc_pre = true;
	%                 	En_err0 = -1;
	%                 	curr_max_it = it + obj.maxIt;
	%             	end
        		end
			end
		end

		if (outerit > obj.OuterLoops || max(abs(err0_array)) < obj.tiny)
			outerStop = true;
		end

	end
    
	obj.physics.Commit("Pathdep");
    obj.physics.Commit("Timedep");
    %figure(9731);
	%semilogy(obj.convergence_log(stepnum,1:it));
	%hold on
	%drawnow();
end

