% Cathy Broderick ILA Feeddback

clear all

%are these changes saved in GIT Desktop? 


min = 6;
duration=100*min*60;
t=1;

% Establish connection with QTM
q =  QMC('QMC_conf_force_10Jan.txt'); %Modify this text to specify the data stream

% Add LabJack
% ljmAsm = NET.addAssembly('LabJack.LJM'); %CB removed labjack 16 Jan 2020
% because it isn't necessary (not sending pulses)
% handle = 0;
% [ljmError, handle] = LabJack.LJM.OpenS('ANY', 'ANY', 'ANY', handle);

% Establish connection with Treadmill and start start treadmill
selfSpeed = input('Input Self-Selected Walking Speed');
setTreadmillSpeed(double(selfSpeed),0.1,double(selfSpeed),'Both');
Paretic_Side = input('Input Participant Paretic (Feedback) Side (R or L)');


Lstp_flag = 0;
Rstp_flag = 0;
% stp_counter = 0;
% memory = 0;
% stp = 1000000;
% ILA = [];

% Pre-Allocate matrices
% GRF_Lz = [];
% GRF_Rz = [];
% R_IC = [];
% R_GTO = [];
% R_LAT_KNEE = [];
% R_LAT_ANKLE = [];
% L_IC = [];
% L_GTO = [];
% L_LAT_KNEE = [];
% L_LAT_ANKLE = [];
% L5_S1 = [];

% begin
while(t<duration);
    % Tell MATLAB what data to get from QTM (via text file)
    [data3D, dataAnalog, dataForce]=QMC(q);      % gets data
    GRF_Lz(t)=dataForce(3,1);     % vertical force of the 1st forceplate
    GRF_Rz(t)=dataForce(3,2);     % vertical force of the 2nd forceplate
    GRF_Raw = [GRF_Lz(t),GRF_Rz(t)];
    
    R_IC(t) = data3D(1,1)./1000;
    R_GTO(t) = data3D(1,2)./1000;
    R_LAT_KNEE(t) = data3D(1,3)./1000;
    R_LAT_ANKLE(t) = data3D(1,4)./1000;
    L_IC(t) = data3D(1,5)./1000;
    L_GTO(t) = data3D(1,6)./1000;
    L_LAT_KNEE(t) = data3D(1,7)./1000;
    L_LAT_ANKLE(t) = data3D(1,8)./1000;
    %     L5_S1(t) = data3D(1,9)./1000;
    
    % Left Swing
    if GRF_Raw(1) <= 80 %if left leg in swing
        if Lstp_flag == 0 %flag at zero means stance
            Lstp_flag = 1; %flag at one means swing
        end
        
        % Left Heel Strike
    elseif GRF_Raw(1) > 80 & GRF_Raw(1) < 2000 %if left leg is in heel strike
        if Lstp_flag == 1
            Lstp_flag = 0;
            
            % normalize to step length = divide by distance between L_LAT_ANKLE and
            % R_LAT_ANKLE %does medial vs lateral matter? Roemmich papers
            % says lateral
            L_STEP_LENGTH = abs(L_LAT_ANKLE-R_LAT_ANKLE); %correct to calculate step length at heel strike?
            
            %left step length why are some positive? and super small?
            %(200mm vs 500mm?)
            
            %Center of Pelvis Calculation: distance between IC markers
            % *** THIS IS AN ABSOLUTE VALUE CALCULATION *** (right?)
            % ** Therefore, center of pelvis is same absolute value
            % calculation for right and left sides...
            CENTER_OF_PELVIS_L = (abs(R_IC - L_IC)./2);
            
            %a = trailing shank
            R_TRAILING_SHANK = (R_LAT_ANKLE(t)-R_LAT_KNEE(t))./(L_STEP_LENGTH);
            
            %b = trailing thigh
            R_TRAILING_THIGH = (R_LAT_KNEE-R_GTO)./(L_STEP_LENGTH);
            
            %c right trailing pelvis
            R_TRAILING_PELVIS = (R_GTO - R_IC)./(L_STEP_LENGTH);
            
            % pelvis rotation
            R_TRAILING_ROTATION = (R_IC - CENTER_OF_PELVIS_L)./(L_STEP_LENGTH);
            
            %g leading thigh
            L_LEADING_THIGH = (L_GTO - L_LAT_KNEE)./(L_STEP_LENGTH); %%changed
            
            %h = leading shank
            L_LEADING_SHANK = (L_LAT_KNEE - L_LAT_ANKLE)./(L_STEP_LENGTH); %%changed
            
            %f
            L_LEADING_PELVIS = (L_IC-L_GTO)./(L_STEP_LENGTH);
            
            %e
            L_LEADING_ROTATION = (CENTER_OF_PELVIS_L - L_IC)./(L_STEP_LENGTH);
            
        end
    end %end left heel strike
    
    % Right Swing
    if GRF_Raw(2) <= 80 %right leg swing
        if Rstp_flag == 0
            Rstp_flag = 1;
        end
        %Right Heel Strike
    elseif GRF_Raw(2) > 80 & GRF_Raw(2) < 2000 %right leg is at heel strike
        if Rstp_flag == 1
            Rstp_flag = 0;
            %compute right step length at heel strike
            R_STEP_LENGTH = abs(R_LAT_ANKLE-L_LAT_ANKLE);
            
            %compute center of pelvis
            CENTER_OF_PELVIS_R = (abs(R_IC - L_IC)./2);
            
            %a = trailing shank
            L_TRAILING_SHANK = (L_LAT_ANKLE-L_LAT_KNEE)./(R_STEP_LENGTH);
            
            %b = trailing thigh
            L_TRAILING_THIGH = (L_LAT_KNEE-L_GTO)./(R_STEP_LENGTH);
            
            %c right trailing pelvis
            L_TRAILING_PELVIS = (L_GTO - L_IC)./(R_STEP_LENGTH);
            
            %d pelvis rotation
            L_TRAILING_ROTATION = (L_IC - CENTER_OF_PELVIS_R)./(R_STEP_LENGTH);
            
            %g leading thigh
            R_LEADING_THIGH = (R_GTO - R_LAT_KNEE)./(R_STEP_LENGTH);
            
            %h = leading shank
            R_LEADING_SHANK = (R_LAT_KNEE - R_LAT_ANKLE)./(R_STEP_LENGTH);
            
            %f
            R_LEADING_PELVIS = (R_IC-R_GTO)./(R_LAT_ANKLE-L_LAT_ANKLE);
            
            %e
            R_LEADING_ROTATION = (CENTER_OF_PELVIS_R - R_IC)./(R_STEP_LENGTH);
            
            % putting this plot here because we are only plotting at RHS
            % try plotting only if R_Leading_Shank exists
            if exist('R_TRAILING_SHANK') & exist('R_TRAILING_THIGH') &... %used to be inside right heel strike
                    exist('R_TRAILING_PELVIS') & ...
                    exist('R_TRAILING_ROTATION') & ...
                    exist('L_LEADING_THIGH') &  ...
                    exist('L_LEADING_SHANK') & ...
                    exist('L_LEADING_PELVIS') & ...
                    exist('L_LEADING_ROTATION') & ...
                    exist('L_TRAILING_SHANK') & ...
                    exist('L_TRAILING_THIGH') & ...
                    exist('L_TRAILING_PELVIS') & ...
                    exist('R_LEADING_SHANK') & ...
                    exist('R_LEADING_PELVIS') & ...
                    exist('R_LEADING_ROTATION') & ...
                    exist('L_TRAILING_ROTATION') & ...
                    exist('R_LEADING_THIGH')
                
                ILA = abs(L_TRAILING_SHANK - R_TRAILING_SHANK) + ...
                    abs(L_TRAILING_THIGH - R_TRAILING_THIGH) ...
                    + abs(L_TRAILING_PELVIS - R_TRAILING_PELVIS) ...
                    + abs(L_TRAILING_ROTATION  - R_TRAILING_ROTATION) ...
                    + abs(L_LEADING_ROTATION - R_LEADING_ROTATION) ...
                    + abs(R_LEADING_PELVIS - L_LEADING_PELVIS) ...
                    + abs(R_LEADING_THIGH - L_LEADING_THIGH) + ...
                    abs(R_LEADING_SHANK - L_LEADING_SHANK); %removed divide by 1000
                
            end
            %             if exist('ILA')
            % designate R vs L side for feedback
            if Paretic_Side == 'R'
                GRF_ant = GRF_Rz
            else
                GRF_ant = GRF_Lz
            end
            %             % is anterior / posterior only Z direction?
            %  designate that feedback is give in the anterior / posterior direction where anterior is (negative?) z direction
                        if GRF_ant < 0
            %             plot(GRF_ant)
            %
            %             GRF,'--ks',...  % end -1 removes extra dot at zero removed (1:end-1)
            %                 'LineStyle','none',...
            %                 'MarkerSize',10,...
            %                 'MarkerEdgeColor','g',...
            %                 'MarkerFaceColor','g');             % plot the data point
            set(gca,'Color','k');
            %                 xlim([length(ILA)-10 length(ILA)+5]); %see if this
            %                 ILA(1:end-1)/1000
            %                 ylim([-0.1 1.5])
            %                 set(gca,'XTick',[],'YTick', []) %can easily remove y and x axis labels
            title('ILA attempt')
            %
            %             end %ends ILA calculation - before plot?
            t=t+1;
            drawnow;
            % goal line at 0
            hold on
            plot(xlim,[0 0], 'LineWidth', 1.5', 'Color', 'r')
            
        end
        
    end
    
end
QMC(q,'quit'); %release the QMC object

