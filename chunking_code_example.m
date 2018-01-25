%     if size(data_proc_rndm,1) > 1000000 % if our table is too large we need to split it up into bitesize chunks
%         
%         chunk_limits = 1:500000:size(data_proc_rndm);
%         number_of_chunks = size(chunk_limits,2);
%         chunk_storage = struct;
%         
%         % split data into manageable chunks
%         for chunkID = 1:number_of_chunks-1
%             chunk_storage.(['Chunk_',num2str(chunkID)]) = data_proc_rndm(chunk_limits(chunkID):chunk_limits(chunkID+1)-1,:);
%         end
%         chunkID = chunkID+1;
%         chunk_storage.(['Chunk_',num2str(chunkID)]) = data_proc_rndm(chunk_limits(chunkID):end,:);
% 
%         % init the output tables
%         ppFFD_rndm_inROI_tmp = [];
%         ppAFD_rndm_inROI_tmp = [];
%         DistancesToN_rndm = [];
%         
%         disp('Pause Point'); % stop here if you want to run this across multiple machines; change ChunkID below as required.
%         
%         % Process each chunk in turn
%         for chunkID = 1:number_of_chunks
%             [chunk_storage.(['ppFFD_',num2str(chunkID)]), chunk_storage.(['ppAFD_',num2str(chunkID)]), chunk_storage.(['DistsToN_',num2str(chunkID)])] = DTF2ParaFunc(chunk_storage.(['Chunk_',num2str(chunkID)]),DTF2vars,QueryXYtemp);
%             ppFFD_rndm_inROI_tmp = vertcat(ppFFD_rndm_inROI_tmp,chunk_storage.(['ppFFD_',num2str(chunkID)]));
%             ppAFD_rndm_inROI_tmp = vertact(ppAFD_rndm_inROI_tmp,chunk_storage.(['ppAFD_',num2str(chunkID)]));
%             DistancesToN_rndm = horzcat(DistancesToN_rndm,chunk_storage.(['DistsToN_',num2str(chunkID)]));
%         end
%         
%         clear chunk_limits number_of_chunks chunkID chunk_storage
%         
%     else
%         % Do it the normal way
%         [ppFFD_rndm_inROI_tmp,ppAFD_rndm_inROI_tmp,DistancesToN_rndm] = DTF2ParaFunc(data_proc_rndm,DTF2vars,QueryXYtemp);
%     end
