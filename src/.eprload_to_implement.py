            case '.spe':
				FileFormat = 'MagnettechBinary'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_MagnettechBinary(FileName)
			
			case '.xml':
				FileFormat = 'MagnettechXML'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_MagnettechXML(FileName)
			
			case '.esr':
				FileFormat = 'ActiveSpectrum'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_ActiveSpectrum(FileName)
			
			case '.dat':
				FileFormat = 'AdaniDAT'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_AdaniDAT(FileName)

			case '.json':
				FileFormat = 'AdaniJSON'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_AdaniJSON(FileName)
			
			case '.eco':
				FileFormat = 'qese/tryscore'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				[Data,Abscissa,Parameters] = eprload_qeseETH(FileName)
			
			case '.plt':
				FileFormat = 'MAGRES'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_MAGRES(FileName)
			
			case '.spk'|'.ref':
				FileFormat = 'VarianETH'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				[Data,Abscissa,Parameters] = eprload_VarianE9ETH(FileName)
			
			case '.d00':
				FileFormat = 'WeizmannETH'
				warn("File type {FileFormat} not yet implemented.")
				if scale:
					warn("Scaling not supported for this filetype.")
				# [Data,Abscissa,Parameters] = eprload_d00WISETH(FileName)
