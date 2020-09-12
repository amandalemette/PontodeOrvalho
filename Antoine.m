classdef (CaseInsensitiveProperties = true) Antoine 
    
    % ANTOINE is an object class for computing pure component saturation
    % conditions using Antoine's equation. An object is created by specifying
    % equation parameters and temperature limits. A piecewise model is
    % formed by specifying equation parameters over disjoint temperature
    % intervals. The object can be created with two different syntaxes
    %
    %   Syntax 1:       S = Antoine(Tmin, Tmax, A, B, C)
    %
    %       Arguments are scalars or vectors. Tmin, Tmax establish
    %       temperature limits. A, B, C are parameters in Antoine's
    %       equation. Attempts to compute results outside of the
    %       temperature limits return NaN.
    %
    %   Syntax 2:       S = Antoine([Tmin, Tmax, A, B, C])
    % 
    %       Like syntax 1 except the parameter are now columns of a matrix.
    %       Each row corresponds to a disjoint temperature interval. These
    %       provides a way to implement a 'piecewise' Antoine equation.
    %
    %   Syntax 3:       S = Antoine.SI(Tmin, Tmax, A, B, C)
    %
    %       For creating objects for data given in SI units. The default
    %       units are [deg C] and [mmHg].
    %
    % OBJECT METHODS
    %
    %   S.Psat(T) 
    %       Computes the saturation pressure using Antoine's equation.
    %       Psat returns NaN for arguments outside of the temperature
    %       limits.
    %
    %   S.Tsat(P) 
    %       Computes the saturate temperature using Antoine's equation.
    %
    %   S.plot(T) 
    %       Prepares a plot of saturation pressure versus temperature. If a
    %       vector T is included, then the saturation pressures are
    %       computed at the indicated temperatures. If omitted, then
    %       saturation pressures are computed for the range S.Tmin to
    %       S.Tmax.
    %
    % STATIC METHODS
    %
    %   Antoine.data
    %       Returns Antoine parameters for a selected set of compounds. The
    %       command
    %
    %          pcd = Antoine.data
    %       
    %       returns a structure with fields corresponding to selected
    %       species.
    %
    %   Antoine.search(str)
    %       Searches the NIST Chemistry Webbook for Antoine parameters
    %       using a search string. Results are displayed in a browser
    %       window.
    %
    %   Antoine.SI(Tmin, Tmax, A, B, C)
    %       Creates an Antoine object from data given in SI units.
    %       Convenient for using Amtoine parameters from the NIST Webbook.
    %
    %   Antoine.tutorial
    %       Displays a tutorial introduction to vapor-liquid equilibrium
    %       calculations using Antoine's equation.
    %
    %
    % USAGE NOTES
    %
    %   1. The assummed units are [deg C] for temperature, and [mmHg] for
    %      pressure.
    %
    %   2. Anonymous functions are can be use to implement a implement
    %      a function compute saturation temperature or pressure. For
    %      example
    %
    %           S = Antoine(0, 60, 8, 1750, 235);
    %           Psat = @(T) S.Psat(T);
    %
    %      creates a function Psat(T) for use in other calculations.
    %
    % EXAMPLES 
    %
    %   1. Plot the saturation pressure of benzene as a function of 
    %      temperature. Using data for benzene (P in mmHg, T in degrees C)
    %
    %          bnz = Antoine(8, 103, 6.90656, 1211.033, 220.790);
    %          bnz.plot
    %
    %   2. Find the normal boiling point of Benzene. 
    %
    %          Tb = fzero(@(T)bnz.Psat(T) - 760, 60, [20,100]);
    %          Tb = 
    %
    %   3. Plot saturation temperatures of water from 0 to 2 atm using a
    %      'piecewise Antoine model'. 
    %
    %          % Retrieve the Antoine parameters for water from the data
    %          % stored within the Antoine class
    %          pcd = Antoine.data;
    %          h2o = pcd.h20;
    %
    %          p = (0:.01:2);
    %          plot(p,h2o.Tsat(760*p);
    %          xlabel('Pressure [atm]');
    %          ylabel('Temperature [C]');
    %          title('Saturation Temperature of Water');
    %
    %   4. Consult the Antoine_tutorial.m for a comprehensive set of
    %      demonstrating the use of the Antoine class for vapor-liquid
    %      equilibrium calculations.
    %
    %  
    % Revision History
    %     December  7, 2010 -- Revised tutorial, added NIST Webbook search
    %                          method, and method to enter data in SI
    %                          units.
    %     December  5, 2010 -- Added tutorial and internal data set.
    %     November 10, 2010 -- Original version uploaded to Matlab Central
    %
    % Author
    %     Jeffrey Kantor
    %     University of Notre Dame

    properties
        Tmin       % Lower temperature limit
        Tmax       % Upper temperature limit
        A          % Parameter A of Antoine's equation
        B          % Parameter B of Antoine's equation
        C          % Parameter D of Antoine's equation
        Pmin       % Lower pressure limit
        Pmax       % Upper pressure limit
    end
    
    methods
        
        % function to create the Antoine object 
        
        function S = Antoine(Tmin, Tmax, A, B, C)
            
            % Interpret input arguments
            
            if nargin == 1
                if size(Tmin,2) ~= 5
                    error('Incorrect arguments to Antoine.');
                end
                [Tmin,Tmax,A,B,C] = deal(Tmin(:,1),Tmin(:,2),Tmin(:,3), ...
                    Tmin(:,4),Tmin(:,5));
            end
            
            if any(Tmin >= Tmax)
                error('Lower temperature limit is greater than higher limit');
            end
            
            % Require disjoint temperature intervals
            
            [Tmin,idx] = sort(Tmin);
            Tmax = Tmax(idx);
            
            for i = 2:length(Tmin)
                if (Tmax(i-1) > Tmin(i))
                    error('Temperature intervals are not disjoint.');
                end
            end
            
            % Data looks good, so put it in the object for later use
            
            A = A(idx);
            B = B(idx);
            C = C(idx);
            
            S.Tmin = Tmin;
            S.Tmax = Tmax;
            S.A = A;
            S.B = B;
            S.C = C;
            
            % Compute lower and upper pressure bounds
            
            S.Pmin = 10.^(A - B./(Tmin + C));
            S.Pmax = 10.^(A - B./(Tmax + C));
            
            % If temperature intervals are continguous, then force the
            % pressure intervals to be continguous as well.
            
            for i = 2:length(S.Pmin)
                if (S.Tmax(i-1)==S.Tmin(i))
                    S.Pmax(i-1) = 0.5*(S.Pmax(i-1)+S.Pmin(i));
                    S.Pmin(i) = S.Pmax(i-1);
                end
            end
            
        end
        
        function P = Psat(obj,T)
            % Compute saturation pressure as a function of temperature
            P = zeros(size(T));
            P(:) = NaN;
            for i = 1:length(obj.Tmin)
                idx = find((T >= obj.Tmin(i)) & (T <= obj.Tmax(i)));
                P(idx) = 10.^(obj.A(i) - obj.B(i)./(T(idx) + obj.C(i)));
            end
        end
      
        function T = Tsat(obj,P)
            % Compute saturation temperature as a function of pressure
            T = zeros(size(P));
            T(:) = NaN;
            for i = 1:length(obj.Pmin)
                idx = find((P >= obj.Pmin(i)) & (P <= obj.Pmax(i)));
                T(idx) = obj.B(i)./(obj.A(i) - log10(P(idx))) - obj.C(i);
            end
        end
        
        function plot(obj,T,vargin)
            % obj.plot(T)  
            % 	  Prepares a plot of saturation pressure versus
            %     temperature. If a vector T is included, then the
            %     saturation pressures are computed at the indicated
            %     temperatures. If omitted, then saturation pressures are
            %     computed for the range S.Tmin to S.Tmax.
            if nargin < 2 || isempty(T)
                Ta = min(obj.Tmin);
                Tb = max(obj.Tmax);
                T = Ta:((Tb-Ta)/200):Tb;
            end
            if nargin < 3
                vargin = '';
            end
            P = obj.Psat(T);
            plot(T,P,vargin);
            xlabel('Temperature [deg C]');
            ylabel('Saturation Pressure [mmHg]');
            title('Saturation Pressure versus Temperature');
        end

    end
    
    methods(Static)
        
        function pcd = data
            % Returns sample data for selected compounds. Temperature units
            % are deg C, pressure units are mmHg.
            
            % Water (Murphy, Table B.4)
            pcd.h2o  = Antoine([  0   60  8.10785 1750.286 235.0; ...
                                 60  150  7.96681 1668.210 228.0]);
                             
            % GASES
            
            % Nitrogen (NIST Webbook)
            pcd.n2   = Antoine(-210.01,-147.15, 6.6113, 264.651, 266.362);            
            % Oxygen (NIST Webbook)
            pcd.o2   = Antoine(-218.8, -118.82, 6.8274, 340.024, 269.006);                        
            % Ammonia (Murphy, Table B.4)
            pcd.nh3  = Antoine( -83,  60, 7.36050, 926.132, 240.17);
                             
            % ALCOHOLS          
            % Methanol (Murphy, Table B.4)
            pcd.meoh = Antoine([-14,  65, 7.89750, 1474.08, 229.13; ...
                                 65, 110, 7.97328, 1515.14, 232.85]);
            % Ethanol (Murphy, Table B.4)
            pcd.etoh = Antoine(  -2, 100, 8.04494, 1554.3, 222.65);
            
            % AROMATICS            
            % Benzene (Murphy, Table B.4)
            pcd.bnz  = Antoine(   8, 113, 6.90656, 1211.033, 220.790);            
            % Toluene (Murphy, Table B.4)
            pcd.tol  = Antoine(   6, 137, 6.95464, 1344.8,   219.48);
            
            % NORMAL ALKANES            
            % Methane (NIST Webbook)
            pcd.CH4  = Antoine.SI(90.99,189.99,3.9895,443.028,-0.49);            
            % Ethane (NIST Webbook)
            pcd.C2H6 = Antoine.SI([ 91.22, 144.13, 4.50706, 791.300,  -6.422; ...
                                   144.13, 199.91, 3.93835, 659.739, -16.719]);            
            % Propane (NIST Webbook)
            pcd.C3H8 = Antoine.SI([ 166.02, 230.5, 4.01158,  834.260, -22.763; ...
                                    230.5,  277.5, 3.98292,  819.296, -24.417; ...
                                    277.5,  360.8, 4.53678, 1149.360,  24.906]);            
            % Butane (NIST Webbook)
            pcd.nC4  = Antoine.SI([ 135.42, 212.89, 4.70812, 1200.475, -13.013; ...
                                    212.89, 272.66, 3.85002,  909.650, -36.146; ...
                                    272.66, 425,    4.35576, 1175.581,  -2.071]);            
            % n-Pentane (Murphy, Table B.4)
            pcd.nC5  = Antoine( -50,  58, 6.85221, 1064.63, 233.01);            
            % n-Hexane (Murphy, Table B.4)
            pcd.nC6  = Antoine( -25,  92, 6.87601, 1171.17, 224.41);            
            % n-Heptane (Murphy, Table B.4)
            pcd.nC7  = Antoine(  -2, 124, 6.89677, 1264.90, 216.54);            
            % n-Octane (Murphy, Table B.4)
            pcd.nC8  = Antoine(  19, 152, 6.91868, 1351.99, 209.15);
            
        end
        
        function search(str)
            % webbook(str)
            %   Searches for Antoine data in the NIST Chemistry Webbook. If
            %   successful, opens a web browswer window showing Antoine
            %   Equation Parameters. 
            %
            %       Antoine.search('ch4');  % Search for Methane
            %
            %   Note that NIST reports Antoine in SI units (K, bar). These
            %   can be converted
            %
            %       Antoine.SI(Tmin,Tmin,A,B,C)
        
            url = ['http://webbook.nist.gov/cgi/cbook.cgi?',...
                   'Name=', strrep(strtrim(str),' ','%20'), ...
                   '&Units=SI&Mask=4&#Thermo-Phase'];
               
            [sout,status] = urlread(url);
            
            if ~status
                error('Failed to read url');
            end
            
            % Get the headline <h1>
            
            h1 = regexp(sout,'<h1>(.*)</h1>','tokens');
            while iscell(h1)
                h1 = h1{1};
            end
            if isempty(h1)
                disp('The page returned was not recognized.');
                return
            end
            
            % If page returned was 'Search Results' then the search was
            % ambiguous. Open the page in a web browser and quit.
            
            if strcmp(h1,'Search Results')
                disp(['Search for ',str,' returned an ambiguous result.']);
                sout = ['text://<html><body>',sout,'</body></html>'];
                web(url);
                return
            end
            
            % Strip out the anchor
            
            h1 = regexp(h1,'<a(.*)>(.*)</a>','tokens');
            h1 = h1{1}{end};
            
            % Strip out the Antoine Equation Parameters
 
            a = regexp(sout,'<h3>Antoine Equation Parameters</h3>');
            sout = sout(a:end);
            b = 8 + regexp(sout,'</table>');
            sout = sout(1:min(b(:)));
            
            if isempty(sout)
                disp(['No Antoine parameters were found for ',str]);
                return
            end
            
            % Strip out the Javascript plotting reference
            
            [s e] = regexp(sout,'</p>\s+<p>(.*)</p>');
            sout = [sout(1:s-1),'</p>',sout(e+1:end)];
            
            % Display Antoine Parameter Data recovered from Webbook
                    
            web(['text://<html><body><h1>',h1,'</h1>', ...
                sout,'</body></html>']);

        end
        
        function obj = SI(Tmin,Tmax,A,B,C)
            % Antoine.SI
            %   A static method to create an Antoine object using SI
            %   parameters. The syntax is the same as Antoine, i.e.,
            % 
            %     obj = Antoine.SI(Tmin,Tmax,A,B,C)
            %     obj = Antoine.SI([Tmin,Tmax,A,B,C])
            %
            %   Temperature is given in degrees Kelvin, pressure in bars.
            %   The units are converted to degrees Celsius and mmHg.
            
            if nargin == 1
                obj = Antoine(Tmin);
            elseif nargin == 5
                obj = Antoine(Tmin,Tmax,A,B,C);
            else
                error('Incorrect number of arguments');
            end

            obj.A = obj.A + log10(760/1.01325);
            
            obj.C = obj.C + 273.15;
            obj.Tmin = obj.Tmin - 273.15;
            obj.Tmax = obj.Tmax - 273.15;
            
            obj.Pmin = obj.Pmin*760/1.01325;
            obj.Pmax = obj.Pmax*760/1.01325;
        end
        
        
        function tutorial
            % Display a published version of the Antoine tutorial
            
            web('html/Antoine_tutorial.html');
        end
        
    end
   
end
        