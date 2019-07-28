classdef helperModClassFrameStore < handle
%helperModClassFrameStore Manage data for modulation classification
%   FS = helperModClassFrameStore creates a frame store object, FS, that
%   stores the complex baseband signals in a format usable in machine
%   learning algortihms.
%   
%   FS = helperModClassFrameStore(MAXFR,SPF,LABELS) creates a frame store
%   object, FH, with the maximum number of frames, MAXFR, samples per
%   frame, SPF, and expected labels, LABELS.
%   
%   Methods:
%   
%   add(FS,FRAMES,LABEL) adds frame(s), FRAMES, with label, LABEL, to the
%   frame store.
%
%   [FRAMES,LABELS] = get(FS) returns stored frames and corresponding
%   labels from frame store, FS.
%   
%   See also ModulationClassificationWithDeepLearningEaxample.

  properties
    OutputFormat = FrameStoreOutputFormat.IQAsRows
  end
  
  properties (SetAccess=private)
    %NumFrames Number of frames in the frame store
    NumFrames = 0
    %MaximumNumFrames Capacity of frame store
    MaximumNumFrames
    %SamplesPerFrame Samples per frame
    SamplesPerFrame
    %Labels Set of expected labels
    Labels
  end
  
  properties (Access=private)
    Frames
    Label
  end
  
  methods
    function obj = helperModClassFrameStore(varargin)
      %helperModClassFrameStore Store complex I/Q frames
      %    FS = helperModClassFrameStore(MAXFR,SPF,LABELS) returns a frame
      %    store object, FS, to store complex I/Q baseband frames of type
      %    LABEL with frame size of SPF. Frame are stored as a
      %    [SPFxNUMFRAMES] array.
      
      inputs = inputParser;
      addRequired(inputs, 'MaximumNumFrames')
      addRequired(inputs, 'SamplesPerFrame')
      addRequired(inputs, 'Labels')
      parse(inputs, varargin{:})
      
      obj.SamplesPerFrame = inputs.Results.SamplesPerFrame;
      obj.MaximumNumFrames = inputs.Results.MaximumNumFrames;
      obj.Labels = inputs.Results.Labels;
      obj.Frames = ...
        zeros(obj.SamplesPerFrame,obj.MaximumNumFrames);
      obj.Label = repmat(obj.Labels(1),obj.MaximumNumFrames,1);
    end
    
    function add(obj,frames,label,varargin)
      %add     Add baseband frames to frame store
      %   add(FS,FRAMES,LABEL) adds frame(s), FRAMES, with label, LABEL, to
      %   frame store FS.
      
      numNewFrames = size(frames,2);
      if (~isscalar(label) && numNewFrames ~= length(label)) ...
          && (size(frames,1) ~= obj.SamplesPerFrame)
        error(message('comm_demos:helperModClassFrameStore:MismatchedInputSize'));
      end
      
      % Add frames
      startIdx = obj.NumFrames+1;
      endIdx = obj.NumFrames+numNewFrames;
      obj.Frames(:,startIdx:endIdx) = frames;
      
      % Add labels types
      if all(ismember(label,obj.Labels))
        obj.Label(startIdx:endIdx,1) = label;
      else
        error(message('comm_demos:helperModClassFrameStore:UnknownLabel',...
          label(~ismember(label,obj.Labels))))
      end

      obj.NumFrames = obj.NumFrames + numNewFrames;
    end
    
    function [frames,labels] = get(obj)
      %get     Return frames and labels
      %   [FRAMES,LABELS]=get(FS) returns the frames and corresponding
      %   labels in the frame store, FS. 
      %   
      %   If OutputFormat is IQAsRows, then FRAMES is an array of size
      %   [2xSPFx1xNUMFRAMES], where the first row is the in-phase
      %   component and the second row is the quadrature component.
      %   
      %   If OutputFormat is IQAsPages, then FRAMES is an array of size
      %   [1xSPFx2xNUMFRAMES], where the first page (3rd dimension) is the
      %   in-phase component and the second page is the quadrature
      %   component.
      
      switch obj.OutputFormat
        case FrameStoreOutputFormat.IQAsRows
          I = real(obj.Frames(:,1:obj.NumFrames));
          Q = imag(obj.Frames(:,1:obj.NumFrames));
          I = permute(I,[3 1 4 2]);
          Q = permute(Q,[3 1 4 2]);
          frames = cat(1,I,Q);
        case FrameStoreOutputFormat.IQAsPages
          I = real(obj.Frames(:,1:obj.NumFrames));
          Q = imag(obj.Frames(:,1:obj.NumFrames));
          I = permute(I,[3 1 4 2]);
          Q = permute(Q,[3 1 4 2]);
          frames = cat(3,I,Q);
      end
      
      labels = obj.Label(1:obj.NumFrames,1);
    end
    
    function [fsTraining,fsValidation,fsTest] = ...
        splitData(obj,splitPercentages)
      %splitData Split data into training, validation and test
      %   [FSTRAIN,FSVALID,FSTEST]=splitData(FS,PER) splits the stored
      %   frames into training, validation, and test groups based on the
      %   percentages, PER. PER is a three-element vector,
      %   [PERTRAIN,PERVALID,PERTEST], which specifies training,
      %   validation, and test percentages. FSTRAIN, FSVALID, and FSTEST
      %   are the frame stores for training, validation, and test frames.
      
      fsTraining = helperModClassFrameStore(...
        ceil(obj.MaximumNumFrames*splitPercentages(1)/100), ...
        obj.SamplesPerFrame, obj.Labels);
      fsValidation = helperModClassFrameStore(...
        ceil(obj.MaximumNumFrames*splitPercentages(2)/100), ...
        obj.SamplesPerFrame, obj.Labels);
      fsTest = helperModClassFrameStore(...
        ceil(obj.MaximumNumFrames*splitPercentages(3)/100), ...
        obj.SamplesPerFrame, obj.Labels);
      
      for modType = 1:length(obj.Labels)
        rawIdx = find(obj.Label == obj.Labels(modType));
        numFrames = length(rawIdx);
        
        % First shuffle the frames
        shuffleIdx = randperm(numFrames);
        frames = obj.Frames(:,rawIdx);
        frames = frames(:,shuffleIdx);
        
        numTrainingFrames = round(numFrames*splitPercentages(1)/100);
        numValidationFrames = round(numFrames*splitPercentages(2)/100);
        numTestFrames = round(numFrames*splitPercentages(3)/100);
        extraFrames = sum([numTrainingFrames,numValidationFrames,numTestFrames]) - numFrames;
        if (extraFrames > 0)
          numTestFrames = numTestFrames - extraFrames;
        end
        
        add(fsTraining, ...
          frames(:,1:numTrainingFrames), ...
          obj.Labels(modType));
        add(fsValidation, ...
          frames(:,numTrainingFrames+(1:numValidationFrames)), ...
          obj.Labels(modType));
        add(fsTest, ...
          frames(:,numTrainingFrames+numValidationFrames+(1:numTestFrames)), ...
          obj.Labels(modType));
      end
      
      % Shuffle new frame stores
      shuffle(fsTraining);
      shuffle(fsValidation);
      shuffle(fsTest);
    end

    function shuffle(obj)
      %shuffle  Shuffle stored frames
      %   shuffle(FS) shuffles the order of stored frames.
      
      shuffleIdx = randperm(obj.NumFrames);
      obj.Frames = obj.Frames(:,shuffleIdx);
      obj.Label = obj.Label(shuffleIdx,1);
    end
  end
end

