classdef FEM < handle
    %FEM Calculated Finite Element Method
    %   Provided constructors and many methods
    
    properties (Access = public)
       
        Assembler
        Facet
        Edge
        Promoted
        
        TriMesh
        % all qnodes
        Qnodes
        
        Num_nodes
        Num_elems
        Num_edges
        
        Solution
        
        Ref_points
    end
    
    properties (Access = private)
        Ref_mesh
        Mesh
    end
    
    methods
        function this = FEM(edge_points, prec, min_area, PML)
            
            if nargin == 3
                
                PML = [];
                
            end
            
            this.Ref_mesh = Mesh([0 0 1 0 0 1]', 0.5);
            [this.Ref_points, ~, ~, ~, ~] = this.Ref_mesh.promote(prec);
            
            this.Assembler = Assembler();
            
            this.Facet.Integrator = Integrator(2, 2*prec);
            [this.Facet.Qnodes, this.Facet.Weights] = this.Facet.Integrator.export();
            [this.Facet.Ref, this.Facet.RefX, this.Facet.RefY] = ...
                this.Assembler.reference2D(this.Ref_points, this.Facet.Qnodes);
        
            this.Edge.Integrator = Integrator(1, 2*prec);
            [this.Edge.Qnodes, this.Edge.Weights] = this.Edge.Integrator.export();
            [this.Edge.Ref, this.Edge.RefX] = ...
                this.Assembler.reference1D(prec, this.Edge.Qnodes);
              
            if size(edge_points,2) ~= 1
            	this.Mesh =  Mesh(edge_points', min_area, PML');
            else
            	this.Mesh =  Mesh(edge_points, min_area, PML);
            end
            [this.Num_nodes, this.Num_elems, this.Num_edges] = ...
                this.Mesh.export();
            
            
            [this.Promoted.nodes, this.Promoted.elems, this.Promoted.edges, this.Promoted.indices, this.Promoted.neighbors] = ...
                this.Mesh.promote(prec);
            
            this.TriMesh = this.Promoted.elems(1:3, :);
            this.Qnodes  = this.Assembler.qnodes2D(this.Promoted.nodes,...
                this.Facet.Qnodes, this.Promoted.elems);
            
        end
        
        function [Mass] = assema(this, Fcn)
            [I, J , U] = ...
                this.Assembler.assema(this.Promoted.nodes, ...
                this.Promoted.elems, this.Facet.Ref, this.Facet.Weights, Fcn);
            Mass = sparse(I, J ,U);
        end
        
        function [Stiff] = assems(this ,Fcn)
            [I, J ,V] = ...
                this.Assembler.assems(this.Promoted.nodes, ...
                this.Promoted.elems, this.Facet.RefX, this.Facet.RefY,...
                this.Facet.Weights, Fcn);
            
            Stiff = sparse(I,J , V);
        end
        
        function [Robin] = assemlbc(this, Fcn, BC)
            [I, J , W] = ...
                this.Assembler.assemlbc(this.Promoted.nodes, ...
                BC, this.Edge.Ref, this.Edge.Weights, Fcn);
            
            Robin = sparse(I, J ,W, size(this.Promoted.nodes, 2), ...
                size(this.Promoted.nodes, 2) );
        end
        
        function [LoadVector] = asseml(this, Fcn)
            LoadVector = this.Assembler.asseml(this.Promoted.nodes,...
                this.Facet.Qnodes, this.Promoted.elems, this.Facet.Ref, ...
                this.Facet.Weights, Fcn);
        end
        
        function [BCLoadVector] = assemrbc(this, Fcn, BC)
            BCLoadVector = this.Assembler.assemrbc(this.Promoted.nodes,...
                this.Edge.Qnodes, BC, this.Edge.Ref, this.Edge.Weights, Fcn);
        end

        function [GradLoadVector] = assemgradl(this, Fcn_X, Fcn_Y) 
            GradLoadVector = this.Assembler.assemgradl(this.Promoted.nodes,...
                this.Promoted.elems, this.Facet.Ref, this.Facet.RefX, this.Facet.RefY,...
                this.Facet.Weights, Fcn_X, Fcn_Y);
        end
        
        function [MassEnergy] = assem_massenergy(this, qFcn, pFcn)
            [MassEnergy] = this.Assembler.assemble_massenergy(this.Promoted.nodes, this.Promoted.elems, this.Facet.Ref, this.Facet.Weights, qFcn, pFcn);
        end

        function [StiffEnergy] = assem_stiffenergy(this, qFcn, pFcn)
            [StiffEnergy] = this.Assembler.assemble_stiffenergy(this.Promoted.nodes, this.Promoted.elems, this.Facet.RefX, this.Facet.RefY, ...
                            this.Facet.Weights, qFcn, pFcn);
        end

        function [MassEnergyGrad] = assem_massenergygrad(this, qFcn, pFcn)
            [MassEnergyGrad] = this.Assembler.assemble_massenergygrad(this.Promoted.nodes, this.Promoted.elems, this.Facet.Ref, ...
                                this.Facet.Weights, qFcn, pFcn);
        end

        function [StiffEnergyGrad] = assem_stiffenergygrad(this, qFcn, pFcn)
            [StiffEnergyGrad] = this.Assembler.assemble_stiffenergygrad(this.Promoted.nodes, this.Promoted.elems, this.Facet.RefX, ...
                                this.Facet.RefY, this.Facet.Weights, qFcn, pFcn);
        end

        function [MassGrad] = assem_massgrad(this, Fcn_r, Fcn_c)
            MassGrad = this.Assembler.assemble_massgrad(this.Promoted.nodes, this.Promoted.elems, this.Facet.Ref, this.Facet.Weights, Fcn_r, Fcn_c);            
        end
        
        function [StiffGrad] = assem_stiffgrad(this, Fcn_r, Fcn_c)
            StiffGrad = this.Assembler.assemble_stiffgrad(this.Promoted.nodes, this.Promoted.elems, this.Facet.RefX, this.Facet.RefY,...
                this.Facet.Weights, Fcn_r, Fcn_c);            
        end

        function [QWeightVector] = assemb_q2itrans(this)
            [QWeightVector] = this.Assembler.assemble_q2itrans(this.Promoted.nodes, this.Promoted.elems, this.Facet.Ref, this.Facet.Weights);
        end

        function [LoadTransQVector] = assem_loadtrans(this, PVector)
            LoadTransQVector = this.Assembler.assemble_loadtrans(this.Promoted.nodes, this.Promoted.elems, this.Facet.Ref, this.Facet.Weights, PVector);
        end
        
        function [TransPVector] = assem_p2qtrans(this, QVector)
            TransPVector = this.Assembler.assemble_p2qtrans(this.Promoted.nodes, this.Promoted.elems, this.Facet.Ref, QVector);
        end
        
    end
    
    methods(Static)
        
        function [Norm] = norm(u, mass)
            [m, n] = size(u);
            
            if m == 1
                Norm = u*mass*u';
            elseif n == 1
                Norm = u'*mass*u;
            else
                error('FEM:Norm', 'Dimension does not match.\n');
            end
            
            
        end
        
        function [Norm] = semi_norm(u, stiff)
            [m, n] = size(u);
            
            if m == 1
                Norm = u*stiff*u';
            elseif n == 1
                Norm = u'*stiff*u;
            else
                error('FEM:SEMI_NORM', 'Dimension does not match.\n');
            end
            
        end
        
        
    end
    
end
